function C = OutletIXAtNearestStream(FD,S,DEM)

% ChiAtNearestStream computes Chi values at the nearest stream.
%
% Syntax
%
%     C = OutletIXAtNearestStream(FD,S,DEM)
%
% Description
%
%     OutletIXAtNearestStream calculates the outlet IX value of each cell in
%     a digital elevation model DEM at the nearest stream cell in S along
%     the flow path in FD.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)s
%
% Output arguments
%
%     C    Outlet IX value at nearest streams (class: GRIDobj)
%

% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016
% Modified by: Kobi Havusha (kobihavu@post.bgu.ac.il)
% Date: 9. September, 2020
% Modified by: Liran Goren (gorenl@bgu.ac.il)
% Date: 7. September, 2022
% 

narginchk(3,3)

validatealignment(S,DEM);
validatealignment(FD,DEM);
C = DEM; % C is an elevation grid.
C.Z = -inf(DEM.size); % C elevation values are -inf.
outlet_ix = streampoi(S,'outlets','ix');
C.Z(outlet_ix) = outlet_ix; % C.Z values at the outlets are outlets.
ix = FD.ix; % givers
ixc = FD.ixc; % receivers
for r = numel(ix):-1:1 % Goes through all the givers (bottom to top).
    C.Z(ix(r)) = max(C.Z(ix(r)),C.Z(ixc(r))); % C values at the givers index is the maximum value between the giver and receiver (bottom to top).
end

