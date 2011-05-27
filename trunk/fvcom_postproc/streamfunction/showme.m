  clf ; clear all


% read nodal data 
  fin     = 'grid.dat' ;
  [vx,vy] = textread(fin,'%f %f');
  m = prod(size(vx));
  coords = zeros(m,2);
  coords(:,1) = vx;
  coords(:,2) = vy;

% read triangle data
  fin     = 'cell.dat' ;
  [i1,i2,i3] = textread(fin,'%d %d %d');
  n = prod(size(i1));
  nv = zeros(n,3);
  nv(:,1) = i1;
  nv(:,2) = i2;
  nv(:,3) = i3;
  

% read streamfunction 
  fin     = ['psi.dat'];
  [psi] = textread(fin,'%f');

% reverse sign if mostly negative
  if(mean(psi) < 0);
     psi = -psi;
  end;


% plot
  patch('Vertices',coords,'Faces',nv,...
           'Cdata',psi,'edgecolor','interp','facecolor','interp');
  colorbar

