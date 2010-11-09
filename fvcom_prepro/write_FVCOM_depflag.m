function write_FVCOM_depflag(depflag,filename,mytitle) 

% Dump spatially-variable deposition flag (depflag) to FVCOM forcing file
%
% function write_FVCOM_depflag(depflag,filename,mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing spatially variable depflag for FVCOM 
%
% INPUT 
%   depflag   = user defined deposition flag (=0, no deposition, =1, deposition) 
%               on the nodes
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute) 
%
% OUTPUT:
%    NetCDF file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_depflag(depflag, 'tst_depflag.nc', 'no deposition in river')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off
subname = 'write_FVCOM_depflag';
fprintf('\n'); fprintf(['begin : ' subname '\n']);

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(~exist('depflag'))
	error('incorrect usage of gen_depflag_file, must provide depflag field')
end;
if(~exist('filename'))
	error('incorrect usage of gen_depflag_file, must provide filename')
end;
if(~exist('title'))
	error('incorrect usage of gen_depflag_file, must provide title field')
end;

% check dimensions
nVerts = prod(size(depflag));
if(nVerts == 0)
	error('dimension of depflag is 0, something is wrong ')
end;

%------------------------------------------------------------------------------
% Dump to depflag NetCDF file
%------------------------------------------------------------------------------
fprintf('Dumping to depflag NetCDF file: \n',filename);
fprintf('Size of depflag array: \n',nVerts);
nc = netcdf(filename,'clobber');
nc.title = mytitle;
nc('node') = prod(size(depflag));
nc{'depflag'}  = ncint('node');
nc{'depflag'}.long_name = 'deposition flag';
nc{'depflag'}.units = '-';
nc{'depflag'}(1:nVerts) = depflag(1:nVerts);
ierr = close(nc);



fprintf(['end   : ' subname '\n'])


