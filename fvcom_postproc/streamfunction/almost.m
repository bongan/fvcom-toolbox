% generate streamfunction
clear all; close all;


segnames = {'Mainland','Naushon','Pasque','Nashawena','Cuttyhunk','Penikese','MV','Nant','OBC'};
kmlfiles = {'none','Naushon_Box.kml','Pasque_Box.kml','Nashawena_Box.kml','Cuttyhunk_Box.kml','Penikese_Box.kml','MV_Box.kml','Nantucket_box.kml','none'};
ids      = [1,2,3,4,5,6,7];
symbols  = {'k+','ro','g+','b+','mo','co','y+'}; 
nsegs = numel(ids); 
for i=1:nsegs
  seg(i).left    = i-1;
  seg(i).segname = char(segnames(i));
  seg(i).kmlfile = char(kmlfiles(i));
  seg(i).id      = ids(i);
  seg(i).symbol  = char(symbols(i));
  if(seg(i).segname(1
end;

  

file = 'whole.nc';  
obcfile = 'scp4.1_obc.dat';
%file = 'box.nc';
%obcfile = 'none';
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
hold on;
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

allbndry = find(isonb==1);
xbndry = x(allbndry);
ybndry = y(allbndry);

% mark island nodes
for i=2:nsegs-1
  [latb,lonb,dumz] = read_kml(seg(i).kmlfile); 
  [xb,yb] = sp_proj('1802','forward',lonb,latb,'m');
  mark = inpolygon(xbndry,ybndry,xb,yb); 
  seg(i).nodes = allbndry(find(mark==1)); 
  isonb(seg(i).nodes) = seg(i).id; 
end;


% read open boundary nodes
if(~strcmp(obcfile,'none'));
[jnk,obcnodes,junk,junk,junk] = textread(obcfile,'%d %d %d %f %f\n','headerlines',1);
xobc = x(obcnodes); 
yobc = y(obcnodes); 
nobc = numel(obcnodes);
isonb(obcnodes) = seg(end).id;   
seg(end).nodes = obcnodes;


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

% set mainland
seg(1).nodes = find(isonb==1);

% plot coasts
figure
for i=1:nsegs
  plot(x(seg(i).nodes),y(seg(i).nodes),seg(i).symbol); hold on;
end;
legend('coast','naushon','pasque','nash','cutty','peni','obc')


% read vorticity, ua,va
vort = nc{'vorticity'}(frame,:);
ua   = nc{'ua'}(frame,:);
va   = nc{'va'}(frame,:);
zeta = nc{'zeta'}(frame,:);

% set Dirichlet conditions for all solid segments
seg(1).dirichlet = 0.0; %convention mainland to 0
UT = TriScatteredInterp(xc,yc,ua');
VT = TriScatteredInterp(xc,yc,va');
DT = TriScatteredInterp(x,y,zeta'+h);
for i=2:nsegs-1
  % find closest points between left and right segments
  clear dmin;
  xleft = x(seg(i-1).nodes);
  yleft = y(seg(i-1).nodes);
  xrght = x(seg(i).nodes);
  yrght = y(seg(i).nodes);
  for j=1:numel(xleft)   
    dist = sqrt( (xrght-xleft(j)).^2 + (yrght-yleft(j)).^2);
    [dmin(j),imin(j)] = min(dist);
  end;
  [dmin_gl,imin_gl] = min(dmin);  
  ileft = imin_gl; 
  irght = imin(imin_gl);
  xleft = xleft(ileft);
  yleft = yleft(ileft);
  xrght = xrght(irght);
  yrght = yrght(irght);
  plot(xleft,yleft,'ro','MarkerSize',5);
  plot(xrght,yrght,'ko','MarkerSize',5);

  % compute flux through the points 
  flux = 0.;
  nflux_seg  = 10;
  xflux = xleft:(xrght-xleft)/nflux_seg:xrght;
  yflux = yleft:(yrght-yleft)/nflux_seg:yrght;
  for j=1:nflux_seg  
    dx = xflux(j+1)-xflux(j); 
    dy = yflux(j+1)-yflux(j); 
    xm = .5*(xflux(j+1)+xflux(j)); 
    ym = .5*(yflux(j+1)+yflux(j)); 
    u  = UT(xm,ym); 
    v  = VT(xm,ym); 
    d  = DT(xm,ym); 
    %fprintf('accum flux %f %f %f %f %f\n',dx,dy,u,v,d); 
    flux = flux + d*(-u*dy + v*dx);
  end;
  fprintf('flux between %s and %s is %f\n',seg(i-1).segname,seg(i).segname,flux);
  seg(i).dirichlet = seg(i-1).dirichlet + flux;
end;
clear UT
clear VT
clear DT
  

  

% generate bcs in fem2d format
dnodes = [];
dirichlet = [];
for i=1:nsegs-1
  dnodes = [dnodes seg(i).nodes'];
  dirichlet = [dirichlet seg(i).dirichlet*ones(numel(seg(i).nodes),1)'];
end;
nnodes = seg(end).nodes;  
neumann = zeros(numel(nnodes),1);


bd={};
bd{1} = [dnodes' 0.*dirichlet'];
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


