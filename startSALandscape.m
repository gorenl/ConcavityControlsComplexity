function [DEM,drca,mapping] = startSALandscape(len,wid,res)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

% function returns:
% DEM is a GRIDobg
% drca is a data structure, a natrix, that represents the OCN. 
% its columns are = donor | receiver | catchment | darainage area | link length
% mapping is a data structure to map row number to donors. Allows rapid
% search


close all;
global t_area movie reality_check

x = 1:res:len; % column
y = 1:res:wid; % rows
[X,Y] = meshgrid(x,y);
[r,c] = size(X);
Z = rand(r,c)*50+1; %random topography
Z(1,:) = 0; %zero out outlets
Z(end,:) = 0;
Z(:,1) = 0;
Z(:,end) = 0;

DEMi = GRIDobj(X,Y,Z); % Generate GRIDobj 
DEM = fillsinks(DEMi);
FD = FLOWobj(DEM,'preprocess','carve'); %carving remove lakes. 
A  = flowacc(FD);
C = A >= 1;
S = STREAMobj(FD,C);

if reality_check
    plot(S)
    axis equal
    [donx,dony] = ind2coord(DEM,FD.ix);
    [recx,recy] = ind2coord(DEM,FD.ixc);
    figure
    hold on
    axis equal
    plot([donx recx]',[dony recy]','g')
end

% In topotoolbox all boundary nodes drain to one another. 
% This is not a desirable propetry for the current calculation. 
% The following code lines set the reciever of boundary nodes to be themselves. 
boundary_node = DEM.Z(FD.ix) == 0;
don = FD.ix;
rec = FD.ixc;
rec(boundary_node) = don(boundary_node);
%adding the original single outlet to the donor-receiver relation
don = [don; len*wid];
rec = [rec; len*wid];

%extract x-y for donot-receiver 
[donxs,donys] = ind2coord(DEM,don);
[recxs,recys] = ind2coord(DEM,rec);

if reality_check
    % reality check, we can still extract correctly the network  after the last
    % manipulation.
    figure
    hold on
    axis equal
    plot([donxs recxs]',[donys recys]','r')
    for i = 1:length(don)
        if (don(i) == rec(i))
            plot(donxs(i),donys(i),'r*')
        end
    end
end

% order donor-reciever-catchment array from the outlet upward (other LEMs
% call it stack).
drc = zeros(length(don),3);
outlets = find(don == rec); %start with the outlet nodes that drain to themselves

drc(1:length(outlets),1:3) = [don(outlets),rec(outlets),rec(outlets)];
filled = length(outlets); 
point = 1;
while point <= length(don)
    ind = find(rec==drc(point,1));
    if ~isempty(ind)
        for i = 1:length(ind)
            if ~ismember(ind(i),outlets)
                filled = filled+ 1;
                drc(filled,1) = don(ind(i));
                drc(filled,2) = rec(ind(i));
                drc(filled,3) = drc(point,3);
            end
        end
    end
    point = point+1;
end

if reality_check
    %reality check, make sure that the donor-receiver relation was not messed up.
    figure;
    hold on
    axis equal
    for i = 1:length(don)
        [xd,yd] = ind2coord(DEM,drc(i,1));
        [xr,yr] = ind2coord(DEM,drc(i,2));
        plot([xd xr],[yd yr],'c')
    end

    [xd,yd] = ind2coord(DEM,drc(:,1));
    scatter(xd,yd,20,drc(:,3)) % make sure that the catchment assignment is correct.
end

%assgin link length and drainage area in pixels
drca = zeros(length(don),5);
drca(:,1:3) = flipud(drc);
drca(:,4) = 1;
[xd,yd] = ind2coord(DEM,drca(:,1));
[xr,yr] = ind2coord(DEM,drca(:,2));
for i = 1:length(don)
    if drca(i,1) ~= drca(i,2)
        if xd(i) ~= xr(i) && yd(i) ~= yr(i) %diagonal link
            drca(i,5) = sqrt(2*res);
        else %non-diagonal link
            drca(i,5) = res;
        end

        ind = find(drca(:,1) == drca(i,2));
        drca(ind,4) = drca(ind,4) + drca(i,4); %update drainage area of receiver
    end
end

mapping = zeros(length(don),1); % array to map from row number to donor
for i = 1:length(don)
    mapping(drca(i,1)) = i;
end

if reality_check
    figure
    axis equal
    hold on
    for i = 1:length(drca(:,1))
        [xd,yd] = ind2coord(DEM,drca(i,1));
        rec_line = mapping(drca(i,2)); %finding the line of the receiver
        [xr,yr] = ind2coord(DEM,drca(rec_line,1));
        plot([xd xr],[yd yr],'b')
        %text(xd,yd,num2str(drca(i,1)))
    end
    [xd,yd] = ind2coord(DEM,drca(:,1));
    scatter(xd,yd,10,drca(:,4),'filled')
    
end























