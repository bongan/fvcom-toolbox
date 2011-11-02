function write_FVCOM_obs_TS(time,zsl,nverts,tsl,ssl,filename,mytitle) 

% Dump observation profile of T/S to netcdf file to initialize stratification in FVCOM 
%
% function write_FVCOM_obs_TS(mjday,zsl,nverts,tsl,ssl,filename,mytitle) 
%
% DESCRIPTION:
%    Generate a NetCDF file containing vertical profile of T/S for FVCOM 
%
% INPUT 
%   jday= modified julian day or initial model time
%   zsl = zcoordinate of observations, positive up 
%   nverts = number of vertices in the mesh**
%   tsl = temperature at level k (C)
%   ssl = salinity at level k (PSU)
%   filename  = filename to dump to
%   mytitle   = global attribute 
%
% OUTPUT:
%    NetCDF file: filename
%
% **in this script the temp/sal profiles are assumed to be constant at each node
%
% EXAMPLE USAGE
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


% check dimensions
ksl = numel(zsl);

if(numel(tsl) ~= ksl)
  error('dimensions of ssl do not match zsl')
end;
if(numel(ssl) ~= ksl)
  error('dimensions of ssl do not match zsl')
end;

%------------------------------------------------------------------------------
% Dump to S/T profile to NetCDF file 
%------------------------------------------------------------------------------
fprintf('Dumping to NetCDF file: \n',filename);
fprintf('Size of T/S array: \n',ksl);
nc = netcdf(filename,'clobber');
nc.title = mytitle;
nc('ksl') = ksl ;
nc('node') = nverts;
nc('time') = 0;

nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units     = 'days since 0.0';
nc{'time'}.time_zone = 'none';

nc{'Itime'} = ncint('time');
nc{'Itime'}.units     = 'days since 0.0';
nc{'Itime'}.time_zone = 'none';

nc{'Itime2'} = ncint('time');
nc{'Itime2'}.units     = 'msec since 00:00:00';
nc{'Itime2'}.time_zone = 'none';

nc{'zsl'}  = ncfloat('ksl');
nc{'zsl'}.long_name = 'standard z levels positive up';
nc{'zsl'}.units = 'm';

nc{'ssl'}  = ncfloat('time','ksl','node');
nc{'ssl'}.long_name = 'observed_salinity_profile'; 
nc{'ssl'}.units = 'PSU';

nc{'tsl'}  = ncfloat('time','ksl','node');
nc{'tsl'}.long_name = 'observed_temperature_profile'; 
nc{'tsl'}.units = 'C';

% write vars
for i=1:numel(time);
nc{'time'}(i) = time;
nc{'Itime'}(i) = floor(time);
nc{'Itime2'}(i) = mod(time,1)*24*3600*1000.;
end;

nc{'zsl'}(1:ksl) = zsl; 

for i=1:numel(time)
for k=1:ksl
nc{'tsl'}(i,k,:) = tsl(k); 
end;
end;

for i=1:numel(time)
for k=1:ksl
nc{'ssl'}(i,k,:) = ssl(k); 
end;
end;

ierr = close(nc);




 
