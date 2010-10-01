function [lon,lat,x,y] = example_my_project(lon,lat,x,y,direction) 

% Sample user-defined projection and inverse projection of (lon,lat) to (x,y) 
% Copy to my_project (not a member of the toolbox) and modify to suite you
%
% [lon,lat,x,y] = function my_project(lon,lat,x,y)
%
% DESCRIPTION:
%    Define projections between geographical and Euclidean coordinates 
%
% INPUT: 
%   lon       = 1D vector containing longitude
%   lat       = 1D vector containing latitude
%   x         = 1D vector containing x-coordinate
%   y         = 1D vector containing y-coordinate    
%   direction = ['forward' ;  'omverse']
%           
% OUTPUT:
%   (lon,lat) or (x,y) depending on choice of forward or reverse projection
%
% EXAMPLE USAGE
%    [lon,lat,x,y] = my_project(lon,lat,x,y,'reverse') 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'my_project';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

ProjectDirection = 'forward';

if(direction == 'forward')
	ProjectDirection = 'forward';
else
	ProjectDirection = 'inverse';
end;



%------------------------------------------------------------------------------
% Perform the projection:  USER DEFINED 
%------------------------------------------------------------------------------

if(ProjectDirection == 'forward')
	fprintf('Projecting from (lon,lat) to (x,y)\n');
	[x,y] = sp_proj('1802','forward',lon,lat,'m');
	
else
	fprintf('Inverse Projecting from (x,y) to (lon,lat)\n')
	[lon,lat] = sp_proj('1802','inverse',x,y,'m');
end;


fprintf(['end   : ' subname '\n'])


