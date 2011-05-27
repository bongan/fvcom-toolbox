  clf ; clear all
  fout = 'psi.dat';

% calculate streamfunction from forcing field, model domain, and bcs
%
% solve d2(psi)/dx2 = -grad X u 
%       with bcs:  psi = 0.       on solid boundary (dirichlet)
%                  d(psi)/dn = 0. on open boundary  (neumann)
%       
% uses fem2.m: Copyright (c) 2002-04-28, B. Rasmus Anthin.
%
% requires:  grid.dat    from dump_vort 
% requires:  cell.dat    from dump_vort 
% requires:  force.dat   from dump_vort
% requires:  dirich.dat from  dump_vort
% requires:  neumann.dat from dump_vort
% output  :  psi(m) --> streamfunction at nodes

time = greg2mjulian('');
fname = 


%----------------------------------------------------------
% read mesh and determine boundary nodes
%----------------------------------------------------------
x  = nc{'x'}(:);
y  = nc{'y'}(:);
nv = nc{'nv'}(:,:)';
ntsn = nc{'ntsn'}(:);

nbndry 
bndry_nodes
%----------------------------------------------------------
% read open boundary and set Neumann nodes 
%----------------------------------------------------------


  
% read dirichlet boundary nodes
  fin     = 'dirich.dat' ;
  [dnodes,dirich] = textread(fin,'%d %f');
  ndb = prod(size(dnodes));

% read neumann boundary nodes
  fin     = ['neumann.dat'];
  [nnodes,neumann] = textread(fin,'%d %f');
  nnb = prod(size(nnodes));

% read forcing data
  fin     = ['force.dat'];
  [force] = textread(fin,'%f');
  fprintf('data read in \n');

%%%dirichlet boundary condition for analytical test case
%  dirich = zeros(ndb,1);
%  xrange = max(mesh.nodexy(:,1)) - min(mesh.nodexy(:,1));
%  yrange = max(mesh.nodexy(:,2)) - min(mesh.nodexy(:,2));
%  for i=1:ndb
%     xloc = mesh.nodexy(dnodes(i),1);
%     yloc = mesh.nodexy(dnodes(i),2);
%     dirich(i) = 1000.*cos(8*pi*xloc/xrange)+1000.*sin(8*pi*yloc/yrange);
%  end;

% generate bcs in fem2d format
  bd={};
  bd{1} = [dnodes dirich];
  if(nnb > 0)
    bd{2} = [nnodes neumann zeros(nnb,1)];
  end;



% set PDE parameter values
  alpha = ones(m,1);
  beta  = zeros(m,1);
  s     = force;  

% solve 
  fprintf('calling Poisson solve \n');
  [psi] = fem2([x,y],nv,alpha,beta,s,bd);

  patch('Vertices',[x,y],'Faces',nv,...
           'Cdata',psi,'edgecolor','interp','facecolor','interp');
  colorbar

% reverse sign so that streamfunction is predominantly > 0
  if(mean(psi) < 0);
     psi = -psi;
  end;

