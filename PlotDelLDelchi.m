load('DAC_Data_Apr_2023.mat')
theta = DAC_DL_Dchi(:,1);
median_delL = DAC_DL_Dchi(:,3);
median_delchi = DAC_DL_Dchi(:,6);
L_err_pos = abs(DAC_DL_Dchi(:,4)-DAC_DL_Dchi(:,3));
L_err_neg = abs(DAC_DL_Dchi(:,3)-DAC_DL_Dchi(:,2));
chi_err_pos = abs(DAC_DL_Dchi(:,7)-DAC_DL_Dchi(:,6));
chi_err_neg = abs(DAC_DL_Dchi(:,6)-DAC_DL_Dchi(:,5));
figure 
hold on;
errorbar (theta,median_delchi,chi_err_neg,chi_err_pos,'sg')

errorbar (theta,median_delL,L_err_neg,L_err_pos,'sb')

xlabel('\theta ','FontSize',20)
ylabel('normalized \Delta L or \Delta\chi','FontSize',20)
legend('normalized \Delta\chi','normalized \Delta L','FontSize',12)
title('DAC Simulations','FontSize',20)
axis([0 1 0 1.2])