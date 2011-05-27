%function [Mstruc] = extract_mooring(fname,lonpt,latpt,tbeg,tend,deltat)

% Extract a mooring structure from an FVCOM output file 
%
% function [Mstruc] = extract_mooring(fname,lon,lat,tbeg,tend,deltat)
%
% DESCRIPTION:
%    Extract a mooring structure from an FVCOM output file 
%
% INPUT
%    fname : FVCOM output file  
%    lonpt : longitude of mooring location
%    latpt : latitude of mooring location
%    tbeg  : begin time (modified Julian day)
%    tend  : end time (modified Julian day)
%    intv  : frame interval 
%    forces: [false/true] extract forces
%
% OUTPUT:
%    Mstruc = maltab structure containing mooring data 
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

fname = '/Users/gcowles/Projects/GENERIC_MODELS/fvcom/tst/output/tst_0001.nc';
xpt   = 283181; 
ypt   = 154765;
tbeg  = 0.;
tend  = greg2mjulian(1858,11,17,11,0,0);   
oname = 'test.nc';
info  = 'mooring test';

%------------------------------------------------------------------------------
% Set constants
%------------------------------------------------------------------------------
g    = 9.81; %gravitational acceleration

%------------------------------------------------------------------------------
% Read mesh and metrics
%------------------------------------------------------------------------------
if(~exist(fname))
  error([' FVCOM output file ' fname ' does not exist'])
end;

nc = netcdf(fname);
x  = nc{'x'}(:); nVerts = numel(x);
y  = nc{'y'}(:);
xc = nc{'xc'}(:); nElems = numel(xc);
yc = nc{'yc'}(:);
h  = nc{'h'}(:);
tri = nc{'nv'}(:,:)';
nbe = nc{'nbe'}(:,:);
art = nc{'art'}(:);
art1 = nc{'art1'}(:);
siglay = nc{'siglay'}(:,:); 
siglev = nc{'siglay'}(:,:); 
[nLay,jnk] = size(siglay); nLev = nLay + 1;

%[lon,lat] = my_project(x,y,'reverse');
%[xpt,ypt] = my_project(lonpt,latpt);


fprintf('   MESH READ IN  \n');
fprintf('nElems:  %d\n',nElems);
fprintf('nVerts:  %d\n',nVerts);
fprintf('nLay:    %d\n',nLay);  
fprintf('nLev:    %d\n\n',nLev);  

%------------------------------------------------------------------------------
% Find triangle containing point 
%------------------------------------------------------------------------------
dist   = sqrt(  (xc-xpt).^2 + (yc-ypt).^2);

cell = 0;
%try nearest
[rmin,imin] = min(dist);
if(isintriangle(x(tri(imin,1:3)),y(tri(imin,1:3)),xpt,ypt));
  cell = imin;
end;

%sort from min distance to max and search along sort
if(cell ==0)
[distsort,ind] = sort(dist,1,'ascend');
for i=1:nElems
   if(isintriangle(x(tri(ind(i),1:3)),y(tri(ind(i),1:3)),xpt,ypt));
     cell = ind(i);
   end;
end;
end;

if(cell==0);
  error('point %f %f is not in the domain',xpt,ypt);
end;
if(min(nbe(:,cell)) == 0)
  error('point is a boundary point: cannot calculate forces');
end;

fprintf('Point is located in cell %d\n',cell);
fprintf('This is an NOT a boundary point, proceeding\n\n');

%get local metrics
nbe = nbe(:,cell);
aw0 = nc{'aw0'}(:,cell);
awx = nc{'awx'}(:,cell);
awy = nc{'awy'}(:,cell);
a1u = nc{'a1u'}(:,cell);   
a2u = nc{'a2u'}(:,cell);  
nds = tri(cell,:);

figure
patch('Vertices',[x,y],'Faces',tri,...
        'Cdata',h,'edgecolor','k','facecolor','interp');
colorbar
hold on;
plot(xpt,ypt,'ko','MarkerSize',10)
plot(x(nds),y(nds),'r+','MarkerSize',6)
plot(xc(cell),yc(cell),'g+','MarkerSize',6)


%------------------------------------------------------------------------------
% Check times 
%------------------------------------------------------------------------------
time = nc{'time'}(:);
fprintf('begin fvcom time: %s\n',mjul2str(time(1))); 
fprintf('end   fvcom time: %s\n',mjul2str(time(end))); 
fprintf('begin mooring time: %s\n',mjul2str(tbeg)); 
fprintf('end   mooring time: %s\n',mjul2str(tend));  
if(tbeg < time(1)) | (tbeg > time(end)) | (tend < time(1)) | (tend > time(end))
  fprintf('mooring time is not contained within output time bounds\n');
end;
[jnk,ibeg] = min(abs(time-tbeg));
[jnk,iend] = min(abs(time-tend));
fprintf('begin interval %d\n',ibeg);
fprintf('end   interval %d\n',iend);


%==============================================================================
% Main Loop over frames 
%==============================================================================
nco = netcdf(oname,'clobber');
nco.info = info;
nco.from_file = fname;
for f=ibeg:iend



end;

% close up the NetCDF files
close(nc);
close(nco);
  


