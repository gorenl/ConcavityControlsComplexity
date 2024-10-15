
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

% Plotting the data used in the manuscript
%Figs 2c, 4, and 5.
load('SA_160423_010203.mat')
Theta_all = theta;
E_all_all = E_all;
Prob_all_all = Prob_all;
Med_delL_all_all = med_delL_all;
med_delChi_all_all = med_delChi_all;
drca_all_all = drca_all;
load('SA_160423_040506.mat')
Theta_all = [Theta_all , theta];
E_all_all = [E_all_all ; E_all];
Prob_all_all = [Prob_all_all ; Prob_all];
Med_delL_all_all = [Med_delL_all_all;med_delL_all];
med_delChi_all_all = [med_delChi_all_all; med_delChi_all];
drca_all_all = cat(3,drca_all_all , drca_all);
load('SA_160423_070809.mat')
Theta_all = [Theta_all , theta];
E_all_all = [E_all_all ; E_all];
Prob_all_all = [Prob_all_all ; Prob_all];
Med_delL_all_all = [Med_delL_all_all;med_delL_all];
med_delChi_all_all = [med_delChi_all_all; med_delChi_all];
drca_all_all = cat(3,drca_all_all , drca_all);
close all

figure
hold on
subplot(6,1,1)
title('Initial conditions')
plotLandscape(DEM,drca_init,mapping,t_area,1);
k = 2;
for i = length(Theta_all):-2:1
    subplot(6,1,k)
    k = k+1;
    plotLandscape(DEM,drca_all_all(:,:,i),mapping,t_area,1);
    title(strcat('\theta = ',num2str(Theta_all(i))))
end




figure
hold on
n = length(Prob_all_all);
for i = 1:length(Theta_all)
    plot((1:n)/1e6,E_all_all(i,:)/E_all_all(i,1));
    text(n/1e6,E_all_all(i,n)/E_all_all(i,1),...
        strcat('\theta =',num2str(Theta_all(i))));
end
ylabel('P/P_{init}(\theta)')
xlabel('Iterations [\times 10^6]')
title('Energy')


figure
hold on
for i = 1:length(Theta_all)
    plot(1:n,Med_delL_all_all(i,:));
    text(n,Med_delL_all_all(i,n),...
        strcat('\theta =',num2str(Theta_all(i))));
end
ylabel('median \Delta L (\theta)')
xlabel('Iterations [\times 10^6]')
title('median \Delta L')

figure
hold on
for i = 1:length(Theta_all)
    plot(1:n,med_delChi_all_all(i,:));
    text(n,med_delChi_all_all(i,n),...
        strcat('\theta =',num2str(Theta_all(i))));
end
ylabel('median \Delta \chi(\theta)')
xlabel('Iterations [\times 10^6]')
title('\Delta \chi for all tested configurations')

figure
hold on
for i = 1:length(Theta_all)
    plot(Theta_all(i),sum(Prob_all_all(i,:)/n)*1e3,'xb');
    
end
ylabel('\theta')
xlabel('Acceptance ratio [\times 10^-3]')
title('Acceptance ratio')


medDL_vec = zeros(length(theta),1);
q25DL_vec = zeros(length(theta),1);
q75DL_vec = zeros(length(theta),1);
medDChi_vec = zeros(length(theta),1);
q25DChi_vec = zeros(length(theta),1);
q75DChi_vec = zeros(length(theta),1);

for i = 1:length(Theta_all)
    [L,chi] = CalcFlowLenAndChi(drca_all_all(:,:,i),Theta_all(i));
    [medianDelL,q025DelL,q075DelL,medianDelChi,q025DelChi,q075DelChi] = ...
        DeltaLDeltaChi(DEM,drca_all_all(:,:,i),mapping,wid,L,chi);
    medDL_vec(i) = medianDelL;
    q25DL_vec(i) = q025DelL;
    q75DL_vec(i) = q075DelL;
    medDChi_vec(i) = medianDelChi;
    q25DChi_vec(i) = q025DelChi;
    q75DChi_vec(i) = q075DelChi;
end

L_err_pos = abs(q75DL_vec-medDL_vec);
L_err_neg = abs(q25DL_vec-medDL_vec);
chi_err_pos = abs(q25DChi_vec-medDChi_vec);
chi_err_neg = abs(q25DChi_vec-medDChi_vec);
figure 
hold on;
errorbar (Theta_all,medDChi_vec,chi_err_neg,chi_err_pos,'sg')
errorbar (Theta_all,medDL_vec,L_err_neg,L_err_pos,'sb')
xlabel('\theta ','FontSize',20)
ylabel('normalized \Delta L or \Delta\chi','FontSize',20)
legend('normalized \Delta\chi','normalized \Delta L','FontSize',12)
title('OCN Analysis','FontSize',20)
axis([0 1 0 1.2])   
