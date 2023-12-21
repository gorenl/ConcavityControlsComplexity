function res_vector = ConcavityAndComplexityBasedOnPartial(name,minmax,xxyy,nplot)


%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%
% res_vector contains statistics over theta, del_L, and del_chi of the
% elongated range.

close all

resloc = find(name(2:end) == 'P');
rangename = name(1:resloc);
resolution = '145';
%numloc = find(name == '1');
%resolution = name(numloc:numloc+2);
%if isempty(numloc)
%    numloc = find(name == '5');
%    resolution = name(numloc:numloc+1);
%end



load(name)
%defining the outlets belonging to the two sides of the mountain range
outlets1_ind = streampoi(ST_Selected1,'outlets','ix');
outlets2_ind = streampoi(ST_Selected2,'outlets','ix');

figure
% here plot slop - area 
a = getnal(ST_Selected,A)*(DEMc.cellsize^2);
g = gradient(ST_Selected,DEMc);
plot(a,g,'oc')
set(gca,'Xscale','log','Yscale','log');
hold on
sa = slopearea(ST_Selected,DEMc,A);
sa.theta
title(strcat(name, ' \theta =',num2str(-sa.theta)))


figure
hold on
imageschs(DEMc,[],'usepermanent',true,'colormap','gray');

plot(ST_Selected1,'m')
plot(ST_Selected2,'r')
axis equal

FD.fastindexing = true;

outletIX = OutletIXAtNearestStream(FD, ST_Selected, DEMc);

%defininfg the divide network
D = DIVIDEobj(FD,ST_Selected,'type','strahler');

%we need to inspect segments of the divide to define which belong to the 
% main divide, separating ST_Selected1 ftom ST_Selected2
nan_ind = find(isnan(D.IX));
nan_ind = nan_ind - 1;
[~,~,pix,qix] = getvalue(D,DEMc);
MainDivide_ind = [];
num_segments = 0;
%keeping only divide segments that are shared between outlets belonging to
%opposing sides of the mountain range
for i = 1:length(nan_ind)% run over number of segments
    local_index = nan_ind(i);
    outlet_a = outletIX.Z(pix(local_index));
    outlet_b = outletIX.Z(qix(local_index));
    if (ismember(outlet_a,outlets1_ind) && ismember(outlet_b,outlets2_ind)) || ...
            (ismember(outlet_a,outlets2_ind) && ismember(outlet_b,outlets1_ind))
        %need to account for this segment as the main divide
        if i > 1
            MainDivide_ind = [MainDivide_ind nan_ind(i-1)+2:local_index];
        else
            MainDivide_ind = [MainDivide_ind 1:local_index];
        end
        num_segments=num_segments+1;
    end
end

% Extreme point from which to map the main divide
[xx,yy] = ind2coord(D,D.ep);
if xxyy == 1 && minmax == 1
    [~,cpoint] = min(xx);
elseif xxyy == 1 && minmax == 2
    [~,cpoint] = max(xx);
elseif xxyy == 2 && minmax == 1
    [~,cpoint] = min(yy);
else
    [~,cpoint] = max(yy);
end

new_dist = dist2node(D,D.ep(cpoint)); %distance of divide nodes relative to a far away junction.
Ddist = new_dist(MainDivide_ind);  %same distance only over main divide pixels
[Dx,Dy] = ind2coord(D,D.IX(MainDivide_ind)); % these should be the x,y of the main divide
Dxy_dist = unique([Dx,Dy,Ddist],'rows');
Dxy_dist_sorted = sortrows(Dxy_dist,3);
Dx_sorted = Dxy_dist_sorted(:,1);
Dy_sorted = Dxy_dist_sorted(:,2);

%now that the main divide is ordered we can define sinuosity of the main
%divide
euclidean_divide_distance = sqrt((Dx_sorted(1)- Dx_sorted(end))^2+...
    (Dy_sorted(1)- Dy_sorted(end))^2);
along_divide_distance = (length(Dx_sorted)-1)*D.cellsize;
divide_sinuosity_index = along_divide_distance/euclidean_divide_distance;

D2 = divorder(D,'topo');
plot(D2,'limit',[1000 inf],'color','w')

plot(Dx_sorted(1),Dy_sorted(1),'bx')
plot(Dx_sorted(end),Dy_sorted(end),'bx')
title(name)


%Calculate the concavity index
[mn_best,mn_low,mn_high]=ConcavityDisorder(ST_Selected,DEMc,A,90);
mn_vec = [mn_best,mn_low,mn_high];


%% this section is to calculate delta length and delta chi across divide
%flow distance raster
flow_dist = flowdistance(FD);
%assign to each pixel the flow dist of its nearest channel to which
%it drains. Used in Delta L calculation
DatStream = FlowLenAtNearestStream(FD,ST_Selected,DEMc,flow_dist);
%chi attribute list
chi = chitransform(ST_Selected,A,'mn',mn_best,'plot',0);
%assign to each pixel the chi value of its nearest channel to which
%it drains. Used in Delta chi calculation
ChiatStream = ChiAtNearestStream(FD,ST_Selected,DEMc,chi);
ZatStream = ZAtNearestStream(FD,ST_Selected,DEMc);
AatStream = AAtNearestStream(FD,ST_Selected,DEMc,A);

% find the apparent steepness index
chi_on_stream = chi;
z_on_stream = DEMc.Z(ST_Selected.IXgrid);
fitt = fittype('a*x',...
    'independent','x');
fitobj = fit(chi_on_stream,z_on_stream-min(z_on_stream),fitt);
apparent_steepness_index = fitobj.a;

figure;
plot(chi_on_stream,z_on_stream,'o');
hold on
plot(0:ceil(max(chi_on_stream)),fitobj.a*(0:ceil(max(chi_on_stream)))+min(z_on_stream),'c');
title(rangename)


%calculate del L and del chi over shrinked divideObj based on 
% threshold distance of 1000   
Dshrink = shrink(D,FD,1000);
[~,~,pix,qix] = getvalue(Dshrink,DEMc);

pix = pix(~isnan(pix));
qix = qix(~isnan(qix));
delta_len = zeros(length(pix),1);
delta_chi = zeros(length(pix),1);
n = 0;

for i = 1:length(pix)
    %flow length on the two sides of the divide
    p_len = DatStream.Z(pix(i));
    q_len = DatStream.Z(qix(i));
    % chi on the two sides of the divide
    p_chi = ChiatStream.Z(pix(i));
    q_chi = ChiatStream.Z(qix(i));

    outlet_p = outletIX.Z(pix(i));
    outlet_q = outletIX.Z(qix(i));
    if outlet_p == outlet_q %find the top common junction for the two nodes.
        [IXp,~] = flowpathextract(FD,pix(i));
        [IXq,~] = flowpathextract(FD,qix(i));
        com_junctions = IXp(ismember(IXp,IXq));
        top_junction = com_junctions(1);
        p_len = p_len - DatStream.Z(top_junction);
        q_len = q_len - DatStream.Z(top_junction);
        p_chi = p_chi - ChiatStream.Z(top_junction);
        q_chi = q_chi - ChiatStream.Z(top_junction);
    end

    p_z = ZatStream.Z(pix(i));
    q_z = ZatStream.Z(qix(i));

    p_A = AatStream.Z(pix(i));
    q_A = AatStream.Z(qix(i));


    if (ismember(outlet_p,outlets1_ind) && ismember(outlet_q,outlets1_ind)) || ...
            (ismember(outlet_p,outlets2_ind) && ismember(outlet_q,outlets2_ind)) %% only for internal junctions or basins of the same base level


        if p_chi ~=0 && q_chi~=0 && p_len ~=0 && q_len~=0 && ...
                ~isinf(p_chi) && ~isinf(q_chi) &&...
                ~isinf(p_len) && ~isinf(q_len)
            n = n+1;
            if p_z>q_z % if channel heads don't have the same elevation
                chi_q_fix = q_chi + (p_z-q_z)/apparent_steepness_index; % make a chi correction based on steepness index
                delta_chi(n) = abs(chi_q_fix - p_chi)/mean([chi_q_fix,p_chi]);

                curr_slope_q = apparent_steepness_index*q_A^(-mn_best); %make length elevation based on steepness index and local area
                flow_len_q_fix = q_len + abs(p_z-q_z)/curr_slope_q; %dx = dz/slope
                delta_len(n) = ...
                    abs(p_len-flow_len_q_fix)/mean([p_len,flow_len_q_fix]);
            elseif  q_z>p_z
                chi_p_fix = p_chi + (q_z-p_z)/apparent_steepness_index;
                delta_chi(n) = abs(chi_p_fix - q_chi)/mean([chi_p_fix,q_chi]);
                curr_slope_p = apparent_steepness_index*p_A^(-mn_best); %make length elevation based on steepness index and local area
                flow_len_p_fix = p_len + abs(p_z-q_z)/curr_slope_p; %dx = dz/slope
                delta_len(n) = ...
                    abs(q_len-flow_len_p_fix)/mean([q_len,flow_len_p_fix]);
            else
                delta_len(n) = abs(p_len - q_len)/mean([p_len,q_len]);
                delta_chi(n) = abs(p_chi - q_chi)/mean([p_chi,q_chi]);
                end

                if isnan(delta_len(n)) || isnan(delta_chi(n))
                    'stop';
                end

        end
    end
end
delta_len = delta_len(1:n);
delta_chi = delta_chi(1:n);

del_len_vec = [mean(delta_len) std(delta_len) max(delta_len)];
del_chi_vec = [mean(delta_chi) std(delta_chi) max(delta_chi)];

del_len_quantile = [quantile(delta_len,0.25) quantile(delta_len,0.5) quantile(delta_len,0.75)];

del_chi_quantile = [quantile(delta_chi,0.25) quantile(delta_chi,0.5) quantile(delta_chi,0.75)];



res_vector = [divide_sinuosity_index,apparent_steepness_index,mn_vec,del_len_quantile,del_chi_quantile,n];

savename = strcat(rangename,resolution);
save(savename,'DEM','DEMc','A','FD','ST','ST_Selected','ST_Selected1', ...
    'ST_Selected2','D','Dx_sorted','Dy_sorted','euclidean_divide_distance',...
    'along_divide_distance','divide_sinuosity_index',...
    'mn_vec','del_len_quantile','del_chi_quantile','n');
