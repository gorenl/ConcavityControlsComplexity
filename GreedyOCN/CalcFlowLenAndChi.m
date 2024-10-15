function [flow_len,chi] = CalcFlowLenAndChi(drca,theta)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

%calculate flow length and chi along the random channel network

nodes=length(drca);
drcalchimap = zeros(nodes,8);
drcalchimap(:,1:5) = drca;
drcalchimap(:,8) = 1:nodes; % will be used to restore to original order

%sort according to area assures that we scan nodes from outlet up
areasort_drcalchimap = sortrows(drcalchimap,4,'descend'); 

%generate local maps for a rapid receiver access
local_map = zeros(nodes,1);
for i = 1:nodes
    don = areasort_drcalchimap(i,1);
    local_map(don) = i; %the row # of where don is located
end

for i = 1:nodes
    if areasort_drcalchimap(i,1) == areasort_drcalchimap(i,2) %if outlet
        areasort_drcalchimap(i,6) = 0; % flow length is zero
        areasort_drcalchimap(i,5) = 0; % chi is zero
    else
        rec_num = areasort_drcalchimap(i,2);
        %rec_row = find(areasort_drcalchimap(1:i-1,1)==rec_num);
        rec_row = local_map(rec_num);
        L_of_rec = areasort_drcalchimap(rec_row,6);
        Chi_of_rec = areasort_drcalchimap(rec_row,7);
        areasort_drcalchimap(i,6) = L_of_rec + areasort_drcalchimap(i,5);
        areasort_drcalchimap(i,7) = Chi_of_rec + areasort_drcalchimap(i,5)/areasort_drcalchimap(i,4)^theta;
    end
end

originalsort_drcalchimap = sortrows(areasort_drcalchimap,8);
flow_len = originalsort_drcalchimap(:,6);
chi = originalsort_drcalchimap(:,7);



