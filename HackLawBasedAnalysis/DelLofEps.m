%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

close all
theta = [0:0.01:0.49,0.51:0.01:1];
L1 = 50000;
xc = L1/200;
L2min = 0.2*L1/3.8; 
L2 = (L2min:100:L1);
NdelL = zeros(size(L2));
%solving for the ratio of ka

h = 2;
ka = 3;
for i = 1:length(NdelL)
    NdelL(i) = abs(L1 - L2(i))/mean([L1,L2(i)]);
end


res_mat=(length(NdelL)),zeros(length(theta));
for i = 1:length(theta)
    for j = 1:length(L2)
        pow1 = (1-h*theta(i)) ;
        %my_fun = @(eps) D^pow1 - (1/(eps^theta(i)))*(L-D)^(1-h*theta(i)) - ...
        %    xc^pow1*(1-(1/(eps^theta(i))));
        %res = fzero(my_fun,1);      
        res = ((L1^pow1 - xc^pow1)/(L2(j)^pow1 - xc^pow1))^(1/theta(i));
        if isinf(res) || res == 0 || res < 10^-1. || res > 10^1.
            res = NaN;
        end
        res_mat(j,i)=res;
       % pause(0.02)
    end
end
[X,Y] = meshgrid(theta,NdelL);
figure(1)
hold on;
surf(X,Y,log10(res_mat),'EdgeColor','none');
%set(gca,'Xscale','log')
colormap('cool');
view(2)
xlabel('\theta','FontSize',20)
ylabel('Normalized \Delta L','FontSize',20)
set(gca,'FontSize',20)
hcb=colorbar;
hcb.Label.String = "log_{10}(k_{a_1}/k_{a_2})";
hcb.FontSize = 16;
hcb.Label.FontSize = 20;

%solving for h2 - h1
h1 = 2;

res_mat=(length(NdelL)),zeros(length(theta));
for i = 1:length(theta)
    for j = 1:length(L2)
        pow1 = (1-h1*theta(i)) ;
        my_fun = @(h2) (L1^pow1 - xc^pow1)/(pow1) - ...
            ((L2(j))^(1-h2*theta(i)) - xc^(1-h2*theta(i)))/(1-h2*theta(i));
        
        res = fzero(my_fun,h1);      
        
        if isinf(res) || res < 1.5 || res > 2 
            res = NaN;
        end
        res_mat(j,i)=res;
       % pause(0.02)
    end
end
figure
[X,Y] = meshgrid(theta,NdelL);
figure(2);
hold on;
surf(X,Y,h1 - res_mat,'EdgeColor','none');
%surf(X,Y,res_mat);
colormap('cool');
view(2)
xlabel('\theta','FontSize',20)
ylabel('Normalized \Delta L','FontSize',20)
set(gca,'FontSize',20)
hcb=colorbar;
hcb.Label.String = "h_1 - h_2";
hcb.FontSize = 16;
hcb.Label.FontSize = 20;

%% add curves that represent natural linear ranges

% Theta_mountains = [0.19,0.27,0.35,0.49,0.26,0.38,0.31,0.3,0.42,0.51,0.09,0.8,0.42,...
%     0.24,0.58,0.77,0.17];
% 
% Median_DL = [0.192977071, 0.254687577,0.526836753,0.803062022,0.186519712,...
%     0.317346752,0.240536183,0.401591688,0.459435597,1.150435805,0.15206477...
%     1.366682053,0.797123671,0.450709671,0.840398192,1.395991445,0.261834085];
% 
% DL_25 = [0.084949857,0.118223689,0.237951815,0.392162815,0.079933198,...
% 0.156279907,0.112569964,0.195087261,0.225951567,0.603329524,0.070298128,...
% 0.939993978,0.38027142,0.226190202,0.425358728,0.977809802,0.125402216];
% 
% DL_75 = [0.342927858,0.439360306,0.845360339,1.230599821,0.377710536,...
% 0.548034608,0.436851211,0.677921325,0.729766011,1.581604779,0.268724859,...
% 1.729003131,1.215644896,0.722874194,1.324899793,1.665825129,0.454555981];
% 
% [fitres,gof] = fit(Theta_mountains',Median_DL','power1')
figure(1)
plot3(theta,2.12*theta.^1.6,500*ones(size(theta)),'k');
figure(2)
plot3(theta,2.12*theta.^1.6,500*ones(size(theta)),'k');
figure(1)
plot3(theta,2.49*theta.^1.23,500*ones(size(theta)),'k-.');
figure(2)
plot3(theta,2.49*theta.^1.23,500*ones(size(theta)),'k-.');
figure(1)
plot3(theta,1.54*theta.^2.14,500*ones(size(theta)),'k-.');
ylim([0 1.8])
figure(2)
plot3(theta,1.54*theta.^2.14,500*ones(size(theta)),'k-.');
ylim([0 1.8])
