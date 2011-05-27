  function streamfunction(filein,fileout)

  clf

% interpolate streamfunction onto a regular grid and display using matlab 
% save as a .mat object
%
% requires:  grid.dat    from dump_vort 
% requires:  cell.dat    from dump_vort 
% requires:  streamfunction scalar file (from streamfunction.m) 


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
  
% read streamfunction file 
  fin     = filein 
  [streamfunction] = textread(fin,'%f');

% interpolate onto a regular grid
  ndx  = 500;
  xmin = min(coords(:,1));
  xmax = max(coords(:,1));
  ymin = min(coords(:,2));
  ymax = max(coords(:,2));
  xvec = xmin:(xmax-xmin)/ndx:xmax;
  yvec = ymin:(ymax-ymin)/ndx:ymax;
  sfreg= griddata(coords(:,1),coords(:,2),streamfunction,xvec,yvec');
  sf.xC=xvec;
  sf.yC=yvec;
  sf.stream =sfreg;

% plot using contours
  contour(xvec,yvec,sfreg)

% save data to mat file
  save sf fileout
