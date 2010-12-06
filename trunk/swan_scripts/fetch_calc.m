%function [fetch_struct] = fetch_calc(fvcom_mesh,fvcom_bath,lon,lat,zeta_min,zeta_max,dzeta,dtheta,dry_thresh) 
close all; clear all;
fvcom_mesh = 'skg4.3_grd.dat';
fvcom_bath = 'skg4.3_dep.dat';
lon = -122.4721646 ;
lat = 48.3372476; 
zeta_min = 0;
zeta_max = 0;
nZeta = 1;
nTheta = 36;
dry_thresh = 0.1;
% Calculate fetch as a function of free surface height (positive up) and orientation 
%
% function [fetch_struct] = fetch_calc(fvcom_mesh,lon,lat,zeta_min,zeta_max,dzeta,dtheta,dry_thresh) 
%
% DESCRIPTION:
%   Calculate fetch as a function of free surface height (positive up) and orientation 
%
% INPUT 
%   fvcom_mesh  = FVCOM 3.x grid file 
%   fvcom_bath  = FVCOM 3.x bathymetry file 
%   lon         = longitude of point
%   lat         = latitude of point
%   zeta_min    = minimum value of free surface height (m)
%   zeta_max    = maximum value of free surface height (m)
%   nZeta       = number of zeta partitions 
%   nTheta      = number of theta partitions 
%   dry_thresh  = threshold for a point to be dry
%
% OUTPUT:
%   fetch_struct = fetch structure containing fetch as a function of angle and free surface 
%                  height.  Angle is  
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%-------------------------------------------------
% load the mesh metrics
%-------------------------------------------------

% read the Mesh from an FVCOM grid file
Mobj = read_fvcom_mesh(fvcom_mesh); 
x = Mobj.x;
y = Mobj.y;

% read the bathymetry from an FVCOM grid file
h = read_fvcom_bath(fvcom_bath); 

%-------------------------------------------------
% set grids in zeta and theta 
%-------------------------------------------------
theta = 0.:2*pi/nTheta:2*pi;
zeta  = zeta_min:(zeta_max-zeta_min)/nZeta:zeta_max;
if(zeta_min==zeta_max)
  zeta = zeta_min;
  nZeta = 1;
else
  nZeta = nZeta+1;
end;

%-------------------------------------------------
% project observation point (lon,lat) => (x,y) 
%-------------------------------------------------
[xobs,yobs] = my_project(lon,lat,'forward'); 


%-------------------------------------------------
% compute elevation (positive up) of observation point   
%-------------------------------------------------
zobs = -griddata(x,y,h,xobs,yobs); 

%-------------------------------------------------
% loop over angles and depths, compute fetch 
%-------------------------------------------------
fetch = zeros(nTheta+1,nZeta);
for i=1:nZeta
  elev = zeta(i);
  pts = find( (h+elev) > dry_thresh);
  npts = numel(pts);
  dx  = x(pts)-xobs; 
  dy  = y(pts)-yobs; 
  ang = atan2(dy,dx);
  arc = sqrt(dx.^2 + dy.^2);
  for n=1:npts
    [junk,ind] = min(abs(ang(n)-theta));
    if(arc(n) > fetch(ind,i));  fetch(ind,i) = arc(n); end;
  end;
end;

%-------------------------------------------------
% save to a structure 
%-------------------------------------------------

%-------------------------------------------------
% option to plot 
%-------------------------------------------------
%[X,Y] = meshgrid(zeta,theta*180/pi) ;
%pcolor(X,Y,fetch);
%shading interp;
%xlabel('direction');
%ylabel('zeta');
%colorbar
figure
scatter(x,y,5,h); hold on;
plot(xobs+fetch(:,1)'.*cos(theta),yobs+fetch(:,1)'.*sin(theta),'r')

