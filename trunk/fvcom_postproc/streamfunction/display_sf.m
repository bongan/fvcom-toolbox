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




% read psi 
  fin     = ['psi.dat'];
  [psi] = textread(fin,'%f');

% reverse sign if mostly negative
  if(mean(psi) < 0);
     psi = -psi;
  end;



% plot open boundary for reference
%  plot(xbndry/gscale,ybndry/gscale,'k-')  

% plot depth using contours
%  conts = [-10,-40]; 
%  [c,h] = tricont(nv,coords,depth,conts,'k-');

% plot using contours
%  conts = [.2,.1,0.-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1.];
  if(max(psi) < 1.1)
    conts = [.1,.2,.3,.4,.5,.6,.7,.8,.9];
  elseif (max(psi) > 1.5) 
    conts = [.2,.4,.8,1.,1.2,1.4,1.6,1.8,2.];
  else
    conts = [.15,.3,.45,.6,.75,.9,1.05,1.20,1.35,1.5];
  end;
  [c,h] = tricont(nv,coords,psi,conts);
  clabel(c,h)
  colorbar
  axis([850,1350,-200,250])

%  print -depsc2 
%  eps2pdf fileout 

%  figure

