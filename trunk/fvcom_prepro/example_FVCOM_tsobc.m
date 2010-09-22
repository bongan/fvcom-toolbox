function example_FVCOM_tsobc()
% example file for dumping a file to force temperature and salinity at the open b.
%
% function example_FVCOM_tsobc()
%
% DESCRIPTION:
%    Setup a sample FVCOM hydrographic open boundary forcing file
%
% INPUT
%   
% OUTPUT:
%    FVCOM hydrographic open boundary file
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

warning off;


subname = 'example_FVCOM_tsobc';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

fvcom_bathy = 'tst_dep.dat';
fvcom_obc   = 'tst_obc.dat';
tsOBCFile = 'tst_tsobc.nc';

%------------------------------------------------------------------------------
% read in the FVCOM open boundary node data (need node numbers and dimension)
%------------------------------------------------------------------------------
fid = fopen(fvcom_obc,'r');
if(fid  < 0)
  error(['file: ' fvcom_obc ' does not exist']);
end;
C = textscan(fid, '%s %s %s %s %d', 1);
nObcs = C{5};
obc_nodes = zeros(nObcs,1);
fprintf('reading obc file\n');
fprintf('# nodes %d\n',nObc);
for i=1:nObc
  C = textscan(fid, '%d %d %d', 1);
  obc_nodes(i) = C{2};
end;

fprintf('obc reading complete\n');

%------------------------------------------------------------------------------
% read in the FVCOM bathymetry data (need bathymetry on open boundary nodes)
%------------------------------------------------------------------------------
fid = fopen(fvcom_bathy,'r');
if(fid  < 0)
  error(['file: ' fvcom_bathy ' does not exist']);
end;
C = textscan(fid, '%s %s %s %d', 1);
Nverts = C{4};
h = zeros(Nverts,1);
fprintf('reading bathymetry file\n');
fprintf('# nodes %d\n',Nverts);
for i=1:Nverts
  C = textscan(fid, '%f %f %f', 1);
  h(i) = C{3};
end;
fprintf('min depth %f max depth %f\n',min(h),max(h));
fprintf('bathymetry reading complete\n');
fclose(fid);

%--------------------------------------------------------------
% set variables for NetCDF file
%--------------------------------------------------------------

% extract bathymetry at open boundary nodes
h_obc = h(obc_nodes);

% time
time = 0:1:31.;
nTimes = prod(size(time));

% set siglev/siglay
nSiglay = 10;
nSiglev = 11;
inc = 1./real(nSiglay);
siglev = 0:-inc:-1;
for i=1:nSiglay
	siglay(i) = mean(siglev(i:i+1));
end;


% initialize temperature/salinity arrays
temp = zeros(nObc,nSiglay,nTimes);
salt = zeros(nObc,nSiglay,nTimes);

% set variable temperature and salinity
for i=1:nTimes
	obc_temp(i) = 18. + 2.*real(i-1)/nTimes;
	obc_salt(i) = 30. - 5.*real(i-1)/nTimes;
end;

%--------------------------------------------------------------
% dump to netcdf file
%--------------------------------------------------------------

% open boundary forcing
nc = netcdf(tsOBCFile, 'clobber');       

nc.type = 'FVCOM RIVER FORCING FILE' ;
nc.title = 'simple open boundary hydrography test';   
nc.type =  'FVCOM TIME SERIES OBC TS FILE'; 
nc.history = 'example_FVCOM_tsobc';

% dimensions
nc('nobc') = nObc; 
nc('Datestrln') = 26; 
nc('time') = 0; 
nc('siglay') = nSiglay;
nc('siglev') = nSiglev;

% variables
nc{'river_names'} = ncchar('rivers', 'namelen');

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

nc{'obc_nodes'} = ncint('nobc');
nc{'obc_nodes'}.long_name     = 'Open Boundary Node Number';
nc{'obc_nodes'}.grid = 'obc_grid';

nc{'obc_h'} = ncfloat('nobc');
nc{'obc_h'}.long_name     = 'ocean boundary depth';
nc{'obc_h'}.units = 'm';
nc{'obc_h'}.grid  = 'obc_grid';

nc{'obc_siglev'} = ncfloat('siglev','nobc');
nc{'obc_siglev'}.long_name     = 'ocean_sigma/general_coordinate';
nc{'obc_siglev'}.grid  = 'obc_grid';

nc{'obc_siglay'} = ncfloat('siglay','nobc');
nc{'obc_siglay'}.long_name     = 'ocean_sigma/general_coordinate';
nc{'obc_siglay'}.grid  = 'obc_grid';

nc{'obc_temp'} = ncfloat('time','siglay','nobc');
nc{'obc_temp'}.long_name     = 'sea_water_temperature';
nc{'obc_temp'}.units     = 'Celsius';
nc{'obc_temp'}.grid  = 'obc_grid';

nc{'obc_salinity'} = ncfloat('time','siglay','nobc');
nc{'obc_salinity'}.long_name     = 'sea_water_salinity';
nc{'obc_salinity'}.units     = 'PSU';
nc{'obc_salinity'}.grid  = 'obc_grid';

nc{'obc_nodes'}(:) = obc_nodes;
nc{'obc_h'}(:) = obc_h;
for i=1:nObc
	nc{'obc_siglev'}(1:nSiglev,i) = siglev;
	nc{'obc_siglay'}(1:nSiglay,i) = siglay;
end
for i=1:nTimes
	nc{'time'}(i) = time(i);
	nc{'Itime'}(i) = floor(time(i));
	nc{'Itime2'}(i) = mod(time(i),1)*24*3600*1000.;
	nc{'obc_temp'}(i,:,:) = obc_temp(i);
	nc{'obc_salinity'}(i,:,:) = obc_salt(i);
end;

nc = close(nc);    


fprintf(['end   : ' subname '\n'])
