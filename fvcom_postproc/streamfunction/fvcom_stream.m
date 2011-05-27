% generate streamfunction
clear all; close all;


%file = 'subdom.nc';
%obcfile = 'scp4.1_obc.dat';
file = 'box.nc';
file = 'box_NB_2_Naushon.nc';
file = '/Volumes/Data/Users/gcowles/Projects/NOAA_WHOISG_Waves/idealized_modeling/bb_box/output/box_parabolic.nc'; 
file = '/Volumes/Data/Users/gcowles/Projects/NOAA_WHOISG_Waves/idealized_modeling/bb_box/output/box_NB_2_Naushon.nc';
obcfile = 'none';
nc = netcdf(file);
frame = 1;

% read mesh
x = nc{'x'}(:);
y = nc{'y'}(:);
xc= nc{'xc'}(:);
yc= nc{'yc'}(:);
coords = [x,y];
h = nc{'h'}(:);
nv = nc{'nv'}(:,:)';

patch('Vertices',coords,'Faces',nv,...
           'Cdata',h,'edgecolor','interp','facecolor','interp');
nbe = nc{'nbe'}(:,:)';

% determine all boundary nodes
isonb = zeros(numel(x),1);
for i=1:numel(xc); 
  if(min(nbe(i,1:3))==0)
    if(nbe(i,1) == 0)
      isonb(nv(i,2)) = 1 ; isonb(nv(i,3)) = 1;
    end 
    if(nbe(i,2) ==0) 
      isonb(nv(i,1)) = 1 ; isonb(nv(i,3)) = 1;
    end 
    if(nbe(i,3) ==0) 
      isonb(nv(i,1)) = 1 ; isonb(nv(i,2)) = 1;
    end 
  end;
end;

solid = find(isonb==1);
plot(x(solid),y(solid),'g+')



% read open boundary nodes
if(~strcmp(obcfile,'none'));
[jnk,obcnodes,junk,junk,junk] = textread(obcfile,'%d %d %d %f %f\n','headerlines',1);
xobc = x(obcnodes); 
yobc = y(obcnodes); 
hold on;
plot(xobc,yobc,'r+');
nobc = numel(obcnodes);
isonb(obcnodes) = 2;

% set normals at obc nodes, edges first
for i=2:nobc-1
  dx(i) = xobc(i+1)-xobc(i-1);
  dy(i) = yobc(i+1)-yobc(i-1);
end;
dx(1) = xobc(2)-xobc(1);
dy(1) = yobc(2)-yobc(1);
dx(nobc) = xobc(nobc)-xobc(nobc-1);
dy(nobc) = yobc(nobc)-yobc(nobc-1);
for i=1:nobc
  nx(i) = -dy(i)/sqrt(dx(i)^2 + dy(i)^2);
  ny(i) =  dx(i)/sqrt(dx(i)^2 + dy(i)^2);
end;
scale = .05*(max(x)-min(x));
quiver(xobc,yobc,nx'*scale,ny'*scale);  
axis equal
end;


% read vorticity, ua,va
vort = nc{'vorticity'}(frame,:);
ua   = nc{'ua'}(frame,:);
va   = nc{'va'}(frame,:);

% generate bcs in fem2d format
dnodes = find(isonb==1);
dirich = zeros(numel(dnodes),1); 
nnodes = find(isonb==2);
neumann = zeros(numel(nnodes),1);


bd={};
bd{1} = [dnodes dirich];
if(numel(nnodes) > 0)
  bd{2} = [nnodes neumann zeros(numel(nnodes),1)];
end;

% set PDE parameter values
alpha = ones(numel(vort),1);
beta  = zeros(numel(vort),1);
s     = -vort;

% solve 
  fprintf('calling Poisson solve \n');
  [psi] = fem2([x,y],nv,alpha,beta,s,bd);

  patch('Vertices',[x,y],'Faces',nv,...
           'Cdata',psi,'edgecolor','interp','facecolor','interp');
  colorbar


% reverse sign so that streamfunction is predominantly > 0
%  if(mean(psi) < 0);
%     psi = -psi;
%  end;

% dump psi into netcdf file
close(nc)
nc = netcdf(file,'w');
nc{'psi'} = ncint('time','node');
nc{'psi'}.long_name = 'Streamfunction';
nc{'psi'}.units     = 'm^3/s';
nc{'psi'}(frame,:) = psi;
close(nc)


