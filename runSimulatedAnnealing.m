function [E_all,Prob_all,drca_all,theta] = runSimulatedAnnealing(len,wid,res,stopT,T0_input,alpha_input)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

% A Wrapper code to tun the simulated annealing algorithm for Optimal Channel Network dynamics. 
% The current wrapper focuses on running the same code for different
% concavity index values. It also assumes a greedy, temperature independent
% optimization.
% Parameters: 
% Len and wid: are the domain lentgh and width. 
% res: spacing between grid points. Number of grid points is len/res and
% wid/res
% stopT: Temperature threshold for stopping the iterations OR number of
% global iterations when using the greedy algorithm with a fixed number of
% iterations.
% T0_input: Initial tempeature.
% alpha_input: rate of temperature reduction. Normally 0.8-0.95.



close all;
global t_area movie reality_check
t_area = 5; % thershold drainage area for plotting the OCN
movie = 0; % toggle to generating animations.
reality_check = 0; %toggle for outputting to the screen for debugging

% Generate the initial random network
% drca is the basic data structure, a matrix, to represent the channel
% network
[DEM,drca_init,mapping] = startSALandscape(len,wid,res);


%Alternativiely, previsly generated imnitial netwrok can be loaded
%assumingb it contains the necessary data structures. 
%load('SimulatedAnnealingInitialConditions.mat') 

if reality_check
    plotLandscape(DEM,drca_init,mapping,t_area,1);
end

alpha = alpha_input;

%theta = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
theta = [0.2  0.4  0.6  0.8 ];



T0_max = max(CalculateEnergy(drca_init,theta(1)),T0_input);


%Preparing data structure for output 
E_all = zeros(length(theta),ceil(len*wid/res));
med_delL_all = zeros(length(theta),ceil(len*wid/res));
med_delChi_all = zeros(length(theta),ceil(len*wid/res));
Prob_all = zeros(length(theta),ceil(len*wid/res));
drca_all = zeros(ceil(len*wid/res),5,length(theta));
n_vec = zeros(1,length(theta));
m_vec = zeros(1,length(theta));
T0_vec = zeros(1,length(theta));
lower_energy_vec = zeros(1,length(theta));

% for each final configuration, calculate its energy for all thetas 
final_energy_mat = zeros(length(theta),length(theta));


% saving the initial landscape
save('SimulatedAnnealingInitialConditions','DEM','drca_init','mapping')

for i = 1:length(theta)
    strcat(' \theta = ',num2str(theta(i)))
    T0 = CalculateEnergy(drca_init,theta(i));
    T0_vec(i) = T0;
    [E_vec,Prob_vec,drca,n,m,lower_energy,median_delL_vec,median_delChi_vec] = ...
        simulatedAnnealingLandscape(len,wid,res,stopT,theta(i),T0_max,alpha,DEM,drca_init,mapping);
    E_all(i,1:n) = E_vec;
    med_delL_all(i,1:n) = median_delL_vec;
    med_delChi_all(i,1:n) = median_delChi_vec;
    Prob_all(i,1:n) = Prob_vec;
    drca_all(:,:,i) = drca;
    n_vec(i) = n;
    m_vec(i) = m;
    lower_energy_vec(i) = lower_energy;
    
    plotLandscape(DEM,drca,mapping,t_area,i+1);
    title(strcat('\theta =',num2str(theta(i))))


    % Calculating the energy of the network assuming other theta values. 
    for j = 1:length(theta)
        E_temp = CalculateEnergy(drca,theta(j));
        final_energy_mat(i,j) = E_temp;
    end

end

%Plotting and postprocessing

figure
hold on
for i = 1:length(theta)
    plot(E_all(i,1:n_vec(i))/E_all(i,1));
    text(n_vec(i),E_all(i,n_vec(i))/E_all(i,1),...
        strcat('\theta =',num2str(theta(i)),', lower energy ratio = ',...
        num2str(lower_energy_vec(i)/(ceil(m_vec(i)*wid*len/res)))));
end
ylabel('E/E_{init}(\theta)')
xlabel('Iterations')
title('Emergy for all tested configurations')

figure
hold on
for i = 1:length(theta)

    plot(med_delL_all(i,1:n_vec(i)));
    text(n_vec(i),med_delL_all(i,n_vec(i)),...
        strcat('\theta =',num2str(theta(i))));
end
ylabel('median \Delta L (\theta)')
xlabel('Iterations')
title('median \Delta L all tested configurations')

figure
hold on
for i = 1:length(theta)
    plot(med_delChi_all(i,1:n_vec(i)));
    text(n_vec(i),med_delChi_all(i,n_vec(i)),...
        strcat('\theta =',num2str(theta(i))));
end
ylabel('median \Delta \chi(\theta)')
xlabel('Iterations')
title('normalized median \Delta \chi for all tested configurations')

figure
hold on
for i = 1:length(theta)
    my_ind = logical(Prob_all(i,1:n));
    plot(E_all(i,my_ind)/E_all(i,1));
    my_my_ind = find(my_ind == 1);
    text(lower_energy_vec(i),E_all(i,my_my_ind(end))/E_all(i,1),...
        strcat('\theta =',num2str(theta(i)),', lower energy ratio = ',...
        num2str(lower_energy_vec(i)/(ceil(m_vec(i)*wid*len/res)))));
end
ylabel('E/E_{init}(\theta)')
xlabel('Iterations')
title('Energy of accepted configurations')

% analyzing del_L and del_chi for final optimized OCNs

medDL_vec = zeros(length(theta),1);
q25DL_vec = zeros(length(theta),1);
q75DL_vec = zeros(length(theta),1);
medDChi_vec = zeros(length(theta),1);
q25DChi_vec = zeros(length(theta),1);
q75DChi_vec = zeros(length(theta),1);

for i = 1:length(theta)
    [L,chi] = CalcFlowLenAndChi(drca_all(:,:,i),theta(i));
    [medianDelL,q025DelL,q075DelL,medianDelChi,q025DelChi,q075DelChi] = ...
        DeltaLDeltaChi(DEM,drca_all(:,:,i),mapping,wid,L,chi);
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
errorbar (theta,medDChi_vec,chi_err_neg,chi_err_pos,'sg')
errorbar (theta,medDL_vec,L_err_neg,L_err_pos,'sb')
xlabel('\theta ','FontSize',20)
ylabel('normalized \Delta L or \Delta\chi','FontSize',20)
legend('normalized \Delta\chi','normalized \Delta L','FontSize',12)
title('OCN Analysis','FontSize',20)
axis([0 1 0 1.2])   

%Can save workspace if needed.
