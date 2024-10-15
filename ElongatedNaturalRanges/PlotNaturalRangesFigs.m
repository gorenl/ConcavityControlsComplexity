%Fig. 2a

load('LinearRanges145.mat');
FieldDataMat = FieldDataMat145;
best_theta = FieldDataMat(:,12);
median_delL = FieldDataMat(:,16);
median_delchi = FieldDataMat(:,19);
L_err_pos = abs(FieldDataMat(:,17)-FieldDataMat(:,16));
L_err_neg = abs(FieldDataMat(:,15)-FieldDataMat(:,16));
theta_err_pos = abs(FieldDataMat(:,12)-FieldDataMat(:,14));
theta_err_neg = abs(FieldDataMat(:,12)-FieldDataMat(:,13));
chi_err_pos = abs(FieldDataMat(:,19)-FieldDataMat(:,20));
chi_err_neg = abs(FieldDataMat(:,19)-FieldDataMat(:,18));
figure 
hold on;

errorbar (best_theta,median_delchi,chi_err_neg,chi_err_pos,...
    theta_err_neg,theta_err_pos,'sg')

errorbar (best_theta,median_delL,L_err_neg,L_err_pos,...
    theta_err_neg,theta_err_pos,'sb')

xlabel('\theta ','FontSize',20)
ylabel('normalized \Delta L or \Delta\chi','FontSize',20)
legend('normalized \Delta\chi','normalized \Delta L','FontSize',12)
title('Linear mountain ranges, A_0 = 145 pix','FontSize',20)
axis([0 1 0 2])

%Fig. 6

%% plot data with AI
load('LinearRanges145WithAI.mat')

FieldDataMat = FieldDataMat145;
best_theta = FieldDataMat(:,12);
median_delL = FieldDataMat(:,16);
median_delchi = FieldDataMat(:,19);
L_err_pos = abs(FieldDataMat(:,17)-FieldDataMat(:,16));
L_err_neg = abs(FieldDataMat(:,15)-FieldDataMat(:,16));
theta_err_pos = abs(FieldDataMat(:,12)-FieldDataMat(:,14));
theta_err_neg = abs(FieldDataMat(:,12)-FieldDataMat(:,13));
chi_err_pos = abs(FieldDataMat(:,19)-FieldDataMat(:,20));
chi_err_neg = abs(FieldDataMat(:,19)-FieldDataMat(:,18));
median_AI = FieldDataMat(:,23)/10000;
AI_err_neg = abs(FieldDataMat(:,23)-FieldDataMat(:,22))/10000;
AI_err_pos = abs(FieldDataMat(:,23)-FieldDataMat(:,24))/10000;

index = [(1:length(median_AI))];
figure
errorbar (median_AI(index),median_delL(index),L_err_neg(index),L_err_pos(index),...
    AI_err_neg(index),AI_err_pos(index),'s','Color',[139/256,70/256,0])
xlabel('AI')
ylabel('\Delta L')
