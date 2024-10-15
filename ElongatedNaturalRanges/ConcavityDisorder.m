function [best_mn, mn_low, mn_high] = ConcavityDisorder(ST_Selected,DEMc,A,percent_per_iteration) 
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

%function to find the best suited concavity for StreamObj ST_Selected based
%on the disorder algorithm presented in Gailleton et al. 2021, which is
%based on Hergarten et al., 2016

%Input: ST_Selected is a topotoolbox StreamObj
%          DEMc GridObj for ST_Selected with elevations
%          A GridObj for ST_Selected ith drainage area
%          percent per iteration is percent # of basins to sample in each
%          boot strap iteraction
%Output: best_mn is the best suited concavity index
%             low_mn is a lower bound  concavity index based on
%             bootstrap iteraction
%             high_mn is an upper bound  concavity index based on
%             bootstrap iteraction


mn_vec = 0.05:0.01:1.2; %vec of concavities to check
%mn_vec = 0.1:0.1:1.2; %vec of concavities to check
[L,nc] = conncomps(ST_Selected);
Bootstrap_iterations =ceil(nc*1.5); %number of iterations as the number of basins * 1.5
R = zeros(Bootstrap_iterations,length(mn_vec));
D = zeros(Bootstrap_iterations,length(mn_vec));
Dstar = zeros(Bootstrap_iterations,length(mn_vec));
CS = STREAMobj2cell(ST_Selected);
for j = 1:Bootstrap_iterations
    rand_vec = rand(1,nc);
    ind = find(rand_vec > 1-percent_per_iteration/100); 
    S_bootstrap = union(CS{ind(1)},CS{ind(2)});
    for i = 3:length(ind)
        S_bootstrap = union(S_bootstrap,CS{ind(i)});
    end
    
    for i = 1:length(mn_vec)
        c = chitransform(S_bootstrap,A,'a0',1,'mn',mn_vec(i));
        chi_z_mat = [c, DEMc.Z(S_bootstrap.IXgrid)];
        chi_z_sorted = sortrows(chi_z_mat,2);
        R(j,i) = sum(abs(diff(chi_z_sorted(:,1))));
        D(j,i) = (R(j,i) - chi_z_sorted(end,1))/chi_z_sorted(end,1);
    end
    Dstar(j,:) = D(j,:)/max(D(j,:)); 
end
R_alldata = zeros(1,length(mn_vec));
D_alldata = zeros(1,length(mn_vec));

for i = 1:length(mn_vec)
    c = chitransform(ST_Selected,A,'a0',1,'mn',mn_vec(i));
    chi_z_mat = [c, DEMc.Z(ST_Selected.IXgrid)];
    chi_z_sorted = sortrows(chi_z_mat,2);
    
    %subplot(4,3,i)
    %plot(chi_z_sorted(:,1),chi_z_sorted(:,2))
    %title(num2str(mn_vec(i)))
    R_alldata(i) = sum(abs(diff(chi_z_sorted(:,1))));
    D_alldata(i) = (R_alldata(i) - chi_z_sorted(end,1))/chi_z_sorted(end,1);
end
Dstar_alldata = D_alldata/max(D_alldata);
figure
subplot(1,3,1);
plot(mn_vec,R);
hold on
plot(mn_vec,R_alldata,'k+');
xlabel('mn');
ylabel('R');
subplot(1,3,2);
plot(mn_vec,Dstar);
hold on
plot(mn_vec,Dstar_alldata,'k+');
xlabel('mn');
ylabel('D*');
Dstar_percentile =  prctile(Dstar,[25 50 70]);
subplot(1,3,3);
plot(mn_vec,Dstar_percentile(1,:),'k-.');
hold on
plot(mn_vec,Dstar_percentile(2,:),'r');
plot(mn_vec,Dstar_percentile(3,:),'k');
plot(mn_vec,Dstar_alldata,'k+');
xlabel('mn');
ylabel('D*');
legend('25th percentile', 'median','75th perecntile','all data')
[minDstar, pos_minDstar]= min(Dstar_alldata);
best_mn = mn_vec(pos_minDstar);
ind = find(Dstar_percentile(1,:)< minDstar);
if isempty(ind) %in case the best solution for all data is below the 25 percentile we use the median for the best mn
    [minDstar, pos_minDstar]= min(Dstar_percentile(2,:));
    best_mn = mn_vec(pos_minDstar);
    ind = find(Dstar_percentile(1,:)< minDstar);
end
mn_low = mn_vec(ind(1));
mn_high = mn_vec(ind(end));
