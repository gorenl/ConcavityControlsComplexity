function res_mat = callDeltaChiDeltaLToShare

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%


folders_to_open = {'Theta0.1ASCII';'Theta0.2ASCII';'Theta0.3ASCII';'Theta0.4ASCII';...
    'Theta0.5ASCII';'Theta0.6ASCII';'Theta0.7ASCII';'Theta0.8ASCII';'Theta0.9ASCII'};
file_num = [100000;100000;100000;100000;100000;100000;100000;100000;100000];
theta = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9];
U = [5e-4;5e-4;5e-4;5e-4;5e-4;5e-4;5e-4;5e-4;5e-4];
K = [3.2e-04; 1.5273e-04;2.5936e-05;4.5620e-06;...
    8.3313e-07;1.5811e-07;3.1152e-08;...
    6.3577e-09; 1.3392e-09];
P = [1;1;1;1;1;1;1;1;1];
n = [1;1;1;1;1;1;1;1;1];

mean_del_l = zeros(size(theta));
std_del_l = zeros(size(theta));
max_del_l = zeros(size(theta));
mean_del_chi = zeros(size(theta));
std_del_chi = zeros(size(theta));
max_del_chi = zeros(size(theta));

res_mat = zeros(length(folders_to_open),6);
for i = 1:length(folders_to_open)
    res_mat(i,:)=DeltaChiDeltaLToShare(folders_to_open{i},file_num(i),...
        theta(i),U(i),K(i),P(i),n(i));

end

figure
subplot(2,2,1)
plot(theta,res_mat(:,2),'*-')
xlabel('\theta')
ylabel('median (\Delta L / mean(L_1,L_2))')
subplot(2,2,2)
plot(theta,res_mat(:,3),'*-')
xlabel('\theta')
ylabel('75% quantile (\Delta L / mean(L_1,L_2))')
subplot(2,2,3)
plot(theta,res_mat(:,5),'*-')
xlabel('\theta')
ylabel('median (\Delta \chi / mean(\chi_1,\chi_2))')
subplot(2,2,4)
plot(theta,res_mat(:,6),'*-')
xlabel('\theta')
ylabel('75% quantile (\Delta \chi / mean(\chi_1,\chi_2))')

