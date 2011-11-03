function write_FVCOM_z0(z0,filename,mytitle) 

% Dump spatially-variable bottom roughness (z0) to FVCOM forcing file
%
% function write_FVCOM_z0(z0,filename,mytitle)
%
% DESCRIPTION:
%    Generate a NetCDF file containing spatially variable z0 for FVCOM 
%
% INPUT 
%   z0        = user defined roughness field (m) 
%               roughness is defined on the elements
%   filename  = filename to dump to
%   mytitle   = title of the case (set as global attribute) 
%
% OUTPUT:
%    NetCDF file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_z0(z0field, 'tst_z0.nc', 'z0 tst domain')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off
subname = 'write_FVCOM_z0';
global ftbverbose;
if(ftbverbose);
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(~exist('z0'))
	error('incorrect usage of gen_z0_file, must provide z0 field')
end;
if(~exist('filename'))
	error('incorrect usage of gen_z0_file, must provide filename')
end;
if(~exist('title'))
	error('incorrect usage of gen_z0_file, must provide title field')
end;

% check dimensions
nElems = prod(size(z0));
if(nElems == 0)
	error('dimension of z0 is 0, something is wrong ')
end;

%------------------------------------------------------------------------------
% Dump to z0 NetCDF file
%------------------------------------------------------------------------------
if(ftbverbose);
  fprintf('Dumping to z0 NetCDF file: \n',filename);
  fprintf('Size of z0 array: \n',nElems);
end;
nc = netcdf(filename,'clobber');
nc.title = mytitle;
nc('nele') = prod(size(z0));
nc{'z0b'}  = ncfloat('nele');
nc{'z0b'}.long_name = 'bottom roughness';
nc{'z0b'}.units = 'm';
nc{'z0b'}(1:nElems) = z0(1:nElems);
ierr = close(nc);



if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;


