function [medianDelL,q025DelL,q075DelL,medianDelChi,q025DelChi,q075DelChi] = DeltaLDeltaChi(DEM,drca,mapping,wid,L,chi)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

%find all channel heads
nodes = length(drca);
ind = find(drca(:,4) == 1);
heads_no_sorted = drca(ind,1);
heads = sort(heads_no_sorted);

%the ids of the neighbours of donor i are:
%  i-1-wid   i-1   i-1+wid
%  i-wid      i    i+wid
%  i+1-wid   i+1   1+1+wid

%The algorithm: for each head go over its 4 neighbors with greater index
%for each neighbor check if head, if true calculate delL and delChi
%Works because heads are sorted.
% we do not consider channel heads along the boundaries 
% we do not consider channel heads with length of 1
% we consider only pairs of heads that drain to independent basins. 

% for debuging
%plotLandscape(DEM,drca,mapping,0,3)
%len = nodes/wid;

%Current code performs the calculaiton only on channel heads that belong to
%different catchments. To do the calculaiton over all channel heads needs to 
%uncomment the alternative if statement in each subsection below

delL = [];
delChi = [];
%while length(heads) > 1
for i = 1:length(heads)
    ch1 = heads(i);
    l_ch1 = L(mapping(ch1));
    catch1 = drca(mapping(ch1),3);
    ch1_boundary = DEM.Z(ch1)==0;
    if (~ch1_boundary && l_ch1 > 1)
        chi_ch1 = chi(mapping(ch1));    
        ch2 = ch1-1+wid;
        if  ch2<=nodes
            ch2_boundary = DEM.Z(ch2)==0;
            catch2 = drca(mapping(ch2),3);
            l_ch2 = L(mapping(ch2));
            if (catch2~=catch1) && (~ch2_boundary) && (l_ch2 > 1) && (drca(mapping(ch2),4)==1) %this is a channel head
            % if (~ch2_boundary) && (l_ch2 > 1) && (drca(mapping(ch2),4)==1) %this is a channel head   
                chi_ch2 = chi(mapping(ch2));
                delL = [delL abs(l_ch1-l_ch2)/mean([l_ch1,l_ch2])];
                delChi = [delChi abs(chi_ch1-chi_ch2)/mean([chi_ch1,chi_ch2])];
            end
        end
        ch2 = ch1+wid;
        if  ch2<=nodes
            ch2_boundary = DEM.Z(ch2)==0;
            catch2 = drca(mapping(ch2),3);
            l_ch2 = L(mapping(ch2));
            if  (catch2~=catch1) &&(~ch2_boundary) && (l_ch2 > 1) &&(drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+wid,heads(2:end))
            %if  (~ch2_boundary) && (l_ch2 > 1) &&(drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+wid,heads(2:end))    
                l_ch2 = L(mapping(ch2));
                chi_ch2 = chi(mapping(ch2));
                delL = [delL abs(l_ch1-l_ch2)/mean([l_ch1,l_ch2])];
                delChi = [delChi abs(chi_ch1-chi_ch2)/mean([chi_ch1,chi_ch2])];
            end
        end
        ch2 = ch1+1;
        if  ch2<=nodes
            ch2_boundary = DEM.Z(ch2)==0;
            catch2 = drca(mapping(ch2),3);
            l_ch2 = L(mapping(ch2));
            if (catch2~=catch1) &&(~ch2_boundary) && (l_ch2 > 1) && (drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+1,heads(2:end))
            %if (~ch2_boundary) && (l_ch2 > 1) && (drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+1,heads(2:end))

                
                chi_ch2 = chi(mapping(ch2));
                delL = [delL abs(l_ch1-l_ch2)/mean([l_ch1,l_ch2])];
                delChi = [delChi abs(chi_ch1-chi_ch2)/mean([chi_ch1,chi_ch2])];
            end
        end
        ch2 = ch1+1+wid;
        if  ch2<=nodes
            ch2_boundary = DEM.Z(ch2)==0;
            catch2 = drca(mapping(ch2),3);
            l_ch2 = L(mapping(ch2));
            if (catch2~=catch1) &&(~ch2_boundary) &&  (l_ch2 > 1) &&(drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+1+wid,heads(2:end))
            %if (~ch2_boundary) &&  (l_ch2 > 1) &&(drca(mapping(ch2),4)==1) %this is a channel head %ismember(ch1+1+wid,heads(2:end))
                l_ch2 = L(mapping(ch2));
                chi_ch2 = chi(mapping(ch2));
                delL = [delL abs(l_ch1-l_ch2)/mean([l_ch1,l_ch2])];
                delChi = [delChi abs(chi_ch1-chi_ch2)/mean([chi_ch1,chi_ch2])];
            end
        end
    end
end
medianDelL = quantile(delL,0.50);
q025DelL = quantile(delL,0.25);
q075DelL = quantile(delL,0.75);
medianDelChi = quantile(delChi,0.50);
q025DelChi = quantile(delChi,0.25);
q075DelChi = quantile(delChi,0.75);
