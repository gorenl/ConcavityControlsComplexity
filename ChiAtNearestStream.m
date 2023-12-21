function C = ChiAtNearestStream(FD,S,DEM,chi)

% ChiAtNearestStream computes Chi values at the nearest stream.
%
% Syntax
%
%     C = ChiAtNearestStream(FD,S,DEM,chi)
%
% Description
%
%     ChiAtNearestStream calculates the chi value of each cell in
%     a digital elevation model DEM at the nearest stream cell in S along
%     the flow path in FD.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)
%     chi   node attribute list (nal) of chi values
%
% Output arguments
%
%     C    Chi value at nearest streams (class: GRIDobj)
%
%    Example (For more info, see the attached tutorial file 'DULAB_experiment_Chi')
%
%     DEM = GRIDobj('Diff_EXP_17hr.tif');
%     DEM.Z(DEM.Z<-9998)=NaN;
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     A = flowacc(FD);
%     Ufast = 0.021; % Fast uplift rate is 0.021 m/hr.
%     Umid = 0.016; % Medium uplift rate is 0.016 m/hr.
%     Uslow = 0.008; % Slow uplift rate is 0.008 m/hr.
%     U_K = DEM;
%     U_K.Z = NaN(DEM.size);
%     columns = U_K.size(2);
%     Third = floor(columns/3); % represents 1/3 of the width of the box.
%     U_K.Z(:,1:Third) = Ufast;
%     for i=Third+1:2*Third
%     U_K.Z(:,i) = ((i-Third)/(Third))*(Umid-Ufast)+Ufast;
%     end
%     for i=2*Third+1:columns
%     U_K.Z(:,i) = ((i-2*Third)/(Third))*(Uslow-Umid)+Umid;
%     end
%     chi = ChiPrimeTransform(ST,A,'mn',0.15,'UoverK', U_K);
%     C = ChiAtNearestStream(FD,ST,DEM,chi);
%     imageschs(DEM,C)
%     hold on
%     plot(ST, 'k')
%     hc = colorbar;
%     hc.Label.String = '\chi [m] at nearest stream';
%     
% See also: FLOWobj, FLOWobj/flowdistance, FLOWobj/mapfromnalGRIDobj, 
%           STREAMobj, FLOWobj/vertdistance2stream
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016
% Modified by: Kobi Havusha (kobihavu@post.bgu.ac.il)
% Date: 9. September, 2020
% 

narginchk(4,4)

validatealignment(S,DEM);
validatealignment(FD,DEM);
C = DEM; % C is an elevation grid.
C.Z = -inf(DEM.size); % C elevation values are -inf.
C.Z(S.IXgrid) = chi; % C elevation values at the stream network are chi.
ix = FD.ix; % givers
ixc = FD.ixc; % receivers
for r = numel(ix):-1:1 % Goes through all the givers (bottom to top).
    C.Z(ix(r)) = max(C.Z(ix(r)),C.Z(ixc(r))); % C values at the givers index is the maximum value between the giver and receiver (bottom to top).
end

