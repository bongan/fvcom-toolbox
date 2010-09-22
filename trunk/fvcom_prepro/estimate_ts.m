function [Mobj] = estimate_ts(Mobj) 

% Estimate time step at each node  
%
% [Mobj] = function estimate_ts(Mobj)  
%
% DESCRIPTION:
%    Calculate barotropic time step 
%
% INPUT
%    Mobj = matlab mesh object
%
% OUTPUT:
%    Mobj = matlab structure containing mesh time step
%
% EXAMPLE USAGE
%    Mobj = estimate_ts(Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'estimate_ts';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Set constants
%------------------------------------------------------------------------------

g    = 9.81; %gravitational acceleration
u    = 3.0;  %u-velocity
zeta = 3.0;  %tide amp

if(~Mobj.have_bath)
	error('can''t estimate the time step without bathymetry')
end;

%------------------------------------------------------------------------------
% Compute the time step estimate
%------------------------------------------------------------------------------
x = Mobj.x;
y = Mobj.y;
h = Mobj.h;
tri = Mobj.tri;
nVerts = Mobj.nVerts;
nElems = Mobj.nElems;

ts = ones(nVerts,1)*1e9;
side = zeros(nVerts,1);
for i=1:nElems
  n1 = tri(i,1);
  n2 = tri(i,2);
  n3 = tri(i,3);
  nds = [n1 n2 n3];
  lside = sqrt( (x(n1)-x(n2))^2 + (y(n1)-y(n2))^2); 
  dpth  = max(h(nds))+zeta;
  dpth  = max(dpth,1.);
  ts(nds) = min(ts(nds),lside/( sqrt(g*dpth) + u));
end;
fprintf('minimum time step: %f seconds\n',min(ts))
Mobj.ts = ts;
Mobj.have_ts = true;

fprintf(['end   : ' subname '\n'])
