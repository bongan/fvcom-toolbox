function [Mesh] = read_fvcom_mesh(varargin)

% Read ascii FVCOM mesh files into Matlab mesh object  
%
% [mesh] = function read_fvcom_mesh(varargin)
%
% DESCRIPTION:
%    Read FVCOM grid file, bathymetry file, and boundary node list file 
%    Store in a matlab mesh object 
%
% INPUT [keyword pairs, all optional]:  
%   'grid'       = fvcom grid file [e.g. tst_grd.dat]
%   'bath'       = fvcom bathymetry file [e.g. tst_dep.dat]
%   'bndry'      = fvcom boundary node list file [e.g. tst_obc.dat]
%   'coordinate' = coordinate system [spherical; cartesian (default)]
%   'version'    = FVCOM major version number [2 ; 3 (default)]
%   'project'    = generate (x,y) coordinates if input coordinates are (lon,lat) 
%                  generate (lon,lat) coordinates if input is (x,y)
%                  [true ; false (default)], see project_userdef.m
%
% OUTPUT:
%    Mesh = matlab structure containing mesh data
%
% EXAMPLE USAGE
%    Mesh = read_fvcom_mesh('grid','gom_grd.dat','bath','gom_dep.dat')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

Mobj
Mobj.Nverts = 
