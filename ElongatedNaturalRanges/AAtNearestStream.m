function C = AAtNearestStream(FD,S,DEM,A)

% Based on Kobi Havusha's ChiAtNearestStream 
%
% Syntax
%
%     C = ZAtNearestStream(FD,S,DEM,A)
%
% Description
%
%     AAtNearestStream calculates the drainage area value of each cell in
%     a digital elevation model DEM at the nearest stream cell in S along
%     the flow path in FD.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)
%     A     drainage area (class: GRIDobj)
%
% Output arguments
%
%     C    Flow length value at nearest streams (class: GRIDobj)
%
%   
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016
% Modified by: Kobi Havusha (kobihavu@post.bgu.ac.il)
% Date: 9. September, 2020
% Modified by: Liran Goren (gorenl@bgu.ac.il)
% Date: 7. February, 2022
% 

narginchk(4,4)

validatealignment(S,DEM);
validatealignment(FD,DEM);
C = A; % C is the area raster.
C.Z(~S.IXgrid) = -inf;
ix = FD.ix; % givers
ixc = FD.ixc; % receivers
for r = numel(ix):-1:1 % Goes through all the givers (bottom to top).
    C.Z(ix(r)) = max(C.Z(ix(r)),C.Z(ixc(r))); % C values at the givers index is the maximum value between the giver and receiver (bottom to top).
end

