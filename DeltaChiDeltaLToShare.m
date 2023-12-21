function res_data = DeltaChiDeltaLToShare(folder_name,file_number,m,U,K,P,n)
%this function uses DAC ASCII output files and compute Delta L and Delta
%Chi across divides
%This is used to check the hypothesis that Delta L is a function of theta

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

%close all
%read DAC files
my_str = strcat(folder_name,'/');
net = load(strcat(my_str,'river_network00',num2str(file_number)));
coord = load(strcat(my_str,'coordinates00',num2str(file_number)));
no_channel = load(strcat(my_str,'no_channel_connection00',num2str(file_number)));
tau_z = load(strcat(my_str,'tau_z00',num2str(file_number)));
catchment = load(strcat(my_str,'catchment00',num2str(file_number)));

figure(1)
hold on

%reality check. plot correctly the drainage network
num_node = length(net);
for i = 1:num_node
    if net(i,1) ~=net(i,2)
        xi = coord(i,1);
        yi = coord(i,2);
        xj = coord(net(i,2),1);
        yj = coord(net(i,2),2);
        plot([xi xj],[yi yj],'b');
    end
end
axis equal

%calculate flow length and chi
data_mat = zeros(num_node,4);
data_mat(:,1) = 1:num_node;
data_mat(:,2) = coord(:,3);
sorted_data_mat = sortrows(data_mat,2); %sort by elevation;
for i = 1:num_node
    me = sorted_data_mat(i,1);
    my_rec = net(me,2);
    if my_rec ~= me
        x_rec = coord(my_rec,1);
        y_rec = coord(my_rec,2);
        x_me = coord(me,1);
        y_me = coord(me,2);
        my_area = tau_z(me,4);
        dist = sqrt((x_rec - x_me)^2+(y_rec - y_me)^2);
        ind = find(sorted_data_mat(:,1)== my_rec);
        sorted_data_mat(i,3) = sorted_data_mat(ind,3) + dist; %commulative dist
        sorted_data_mat(i,4) = sorted_data_mat(ind,4) + dist/my_area^m;
    end
end
resort_data_mat = sortrows(sorted_data_mat,1);
flow_len = resort_data_mat(:,3);
chi = resort_data_mat(:,4);

%reality check: flow length and chi increase upstream
%figure
%scatter(coord(:,1),coord(:,2),10,flow_len)
%axis equal
%figure
%scatter(coord(:,1),coord(:,2),10,chi)
%axis equal

delta_chi = [];
delta_l = [];
delta_z= [];
delta_chi_fix = [];
delta_l_fix = [];

%can be used if we want to chek only across the main divide
maxy = max(coord(:,2));
miny = min(coord(:,2));
%steepness index is needed to correct chi and L when the channels 
% heads are at different elevation
steepness_index = (U/(K*P^m))^(1/n);
for i = 1:length(no_channel)
    nodei = no_channel(i,1);
    nodej = no_channel(i,2);
    chi_i = chi(nodei);
    chi_j = chi(nodej);
    if chi_i~=0 && chi_j~=0
        xi = coord(nodei,1);
        yi = coord(nodei,2);
        zi = coord(nodei,3);
        xj = coord(nodej,1);
        yj = coord(nodej,2);
        zj = coord(nodej,3);
        %plot([xi xj],[yi yj],'r');

        catch_i = catchment(nodei);
        y_i = coord(catch_i,2);
        catch_j = catchment(nodej);
        y_j = coord(catch_j,2);
        %if y_i == maxy && y_j == miny || y_j == maxy && y_i == miny
        flow_len_i = flow_len(nodei);
        flow_len_j = flow_len(nodej);
        %if belonging to the same catchment need to find the common junction
        % and remove its flow length and chi value. 
        if catch_i == catch_j 
            node_down_i = [];
            curr_node = nodei;
            while curr_node ~= catch_i
                curr_node = net(curr_node,2);
                node_down_i = [node_down_i curr_node];
            end
            node_down_j = [];
            curr_node = nodej;
            while curr_node ~= catch_j
                curr_node = net(curr_node,2);
                node_down_j = [node_down_j curr_node];
            end
            com_junctions = node_down_i(ismember(node_down_i,node_down_j));
            top_junction = com_junctions(1);
            flow_len_junc = flow_len(top_junction);
            chi_junc = chi(top_junction);
            chi_i =chi_i - chi_junc;
            chi_j =chi_j - chi_junc;
            flow_len_i = flow_len_i - flow_len_junc;
            flow_len_j = flow_len_j - flow_len_junc;
        end
        curr_delta_chi = abs(chi_j-chi_i)/mean([chi_j,chi_i]);
        delta_chi = [delta_chi curr_delta_chi];    
        curr_delta_len = abs(flow_len_i-flow_len_j)/mean([flow_len_i,flow_len_j]);
        delta_l = [delta_l curr_delta_len];
                    

        % when channel heads are not exactly at the same height, we perform a
        % chi fix based on the steepness index and a length fix based on
        % slope area relation for steady-state
        if zi>zj
            chi_j_fix = chi_j + (zi-zj)/steepness_index;
            curr_delta_chi_fix = abs(chi_j_fix-chi_i)/mean([chi_j_fix,chi_i]);
            curr_slope_j = steepness_index*tau_z(no_channel(i,2),4)^(-m/n);
            flow_len_j_fix = flow_len_j + abs(zi-zj)/curr_slope_j; %dx = dz/slope
            curr_delta_l_fix = ...
                abs(flow_len_i-flow_len_j_fix)/mean([flow_len_i,flow_len_j_fix]);
        else 
            chi_i_fix = chi_i + (zj-zi)/steepness_index;
            curr_delta_chi_fix = abs(chi_i_fix-chi_j)/mean([chi_i_fix,chi_j]);
            curr_slope_i = steepness_index*tau_z(no_channel(i,1),4)^(-m/n);
            flow_len_i_fix = flow_len_i + abs(zi-zj)/curr_slope_i; %dx = dz/slope
            curr_delta_l_fix = ...
                abs(flow_len_i_fix-flow_len_j)/mean([flow_len_i_fix,flow_len_j]);
        end     
        delta_chi_fix = [delta_chi_fix curr_delta_chi_fix];
        delta_l_fix = [delta_l_fix curr_delta_l_fix];
        curr_delta_z = abs(zi-zj)/mean([zi,zj]);
        delta_z = [delta_z curr_delta_z];
    end
end

%figure
%hist(delta_l,100)
%title('\Delta L')
figure
hist(delta_l_fix,100)
title('\Delta L fix')
figure
%hist(delta_chi,100)
%title('\Delta \chi')
figure
hist(delta_chi_fix,100)
title('\Delta \chi fix')
%figure
%hist(delta_z,100)
%title('\Delta z')

mmean_len = mean(delta_l);
sstd_len = std(delta_l);
mmax_len = max(delta_l);

mmean_len_fix = mean(delta_l_fix);
sstd_len_fix = std(delta_l_fix);
mmax_len_fix = max(delta_l_fix);

mmean_chi = mean(delta_chi);
sstd_chi = std(delta_chi);
mmax_chi = max(delta_chi);

mmean_chi_fix = mean(delta_chi_fix);
sstd_chi_fix = std(delta_chi_fix);
mmax_chi_fix = max(delta_chi_fix);


mmean_z = mean(delta_z);
sstd_z = std(delta_z);
mmax_z = max(delta_z);

%for the hysteresis - it seems that the fix is exentuating the diff between
%different thetas
%res_data = [quantile(delta_l,0.25) quantile(delta_l,0.5) ...
%   quantile(delta_l,0.75)...
%   quantile(delta_chi,0.25) quantile(delta_chi,0.5) ...
%   quantile(delta_chi,0.75)];

%for normal simulaitons
res_data = [quantile(delta_l_fix,0.25) quantile(delta_l_fix,0.5) ...
    quantile(delta_l_fix,0.75)...
    quantile(delta_chi_fix,0.25) ...
    quantile(delta_chi_fix,0.5) ...
    quantile(delta_chi_fix,0.75)];
