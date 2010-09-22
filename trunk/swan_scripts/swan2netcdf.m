function swan2netcdf(matfile,ncfile,basename,first_time,last_time,increment);
% Convert a SWAN output file (Matlab) into a NetCDF file 
%
% function swan2netcdf(matfile,ncfile,basename,first_time,last_time,increment);
%
% DESCRIPTION:
%    read output from unstructured SWAN model (currently 40.82) and
%    dump to a NetCDF file which is far more useful than a Matlab file.
%
% INPUT 
%   matfile  = Unstructured SWAN Matlab file
%   ncfile   = NetCDF file for output
%   basename = prefix for SWAN mesh, bathymetry, connectivity files
%   first_time:  first time frame in Matlab object frame time
%   last_time:   last time frame in Matlab object frame time
%   increment:   increment in seconds
%
% OUTPUT:
%    NetCDF file containing:
%      a.) time in modified Julian day
%      b.) significant wave height (hs)
%      c.) wave direction
%      d.) mesh
%      e.) peak period
%      f.) wave direction
%      g.) U10
%      h.) V10
%
% EXAMPLE USAGE
%   swan2netcdf('gom1.mat','gom1.nc','gom1','20070101_000000','20070131_000000',3600)
%     this converts gom1.mat to gom1.nc using SWAN grid files gom1.ele, gom1.bot
%     and gom1.node from Jan 1, 2007 00:00:00 to Jan 31, 2007 00:00:00 in increments
%     of 1 hour.
%
% NOTE
%    routine is not refined, e.g. will not check if files exist and will
%    probably crash if you do not have the variables above in the SWAN
%    output file.  
%    You will need approximately the following BLOCK command in your SWAN runfile
%
%    BLOCK 'COMPGRID' NOHEAD 'gom1.mat' LAY 3 XP YP DEP HS RTP TPS DIR WLEN &
%                WINDV OUTPUT 20070101_000000 3600 SEC
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

eval(['load ' matfile]);              % load binary file containing SWAN results
                                      % obtained using BLOCK command with COMPGRID-set
% load the connectivity
elefile=[basename '.ele'];
fid = fopen(elefile);                 % load TRIANGLE element based connectivity file
[nele] = fscanf(fid,'%i',[1 3]);     % get number of triangles
jnk = fscanf(fid,'%i',[4 nele(1)])'; % get connectivity table
tri = jnk(:,2:4);

% load the bathymetry
bathfile=[basename '.bot'];
h = textread(bathfile,'%f\n');
node = prod(size(h));

% check if variable range exists
f = first_time;
l = last_time;
tbeg = datenum(str2num(f(1:4)),str2num(f(5:6)),str2num(f(7:8)),str2num(f(10:11)),str2num(f(12:13)),str2num(f(14:15)));
tend = datenum(str2num(l(1:4)),str2num(l(5:6)),str2num(l(7:8)),str2num(l(10:11)),str2num(l(12:13)),str2num(l(14:15)));
tinc  = increment/(24*3600);
times = tbeg:tinc:tend;
ntimes = prod(size(times));
for i=1:ntimes; 
%  datestr(times(i))
  date = datestr(times(i),30);   
  vname = ['Hsig_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0) 
    error('variable frame %s\n does not exist',vname)
  end;
end;


% dump the netcdf file
nc = netcdf(ncfile, 'clobber');


% dimensions
nc('nele') = nele;
nc('node') = node;
nc('three') = 3;
nc('time') = 0;

% variables

nc{'x'} = ncfloat('node');
nc{'x'}.long_name = 'nodal x-coordinate';
nc{'x'}.units     = 'm'; 

nc{'y'} = ncfloat('node');
nc{'y'}.long_name = 'nodal y-coordinate';
nc{'y'}.units     = 'm'; 

nc{'h'} = ncfloat('node');
nc{'h'}.long_name = 'Bathymetry';  
nc{'h'}.units     = 'm'; 

nc{'d'} = ncfloat('time','node');
nc{'d'}.long_name = 'Depth';  
nc{'d'}.units     = 'm'; 

nc{'nv'} = ncint('three','nele');
nc{'nv'}.long_name = 'nodes surrounding element';

nc{'time'} = ncfloat('time');
nc{'time'}.long_name = 'Modified Julian Day';    
nc{'time'}.units     = 'days'; 

nc{'hs'} = ncfloat('time','node');
nc{'hs'}.long_name = 'Significant Wave Height';
nc{'hs'}.units     = 'm'; 

nc{'wdir'} = ncfloat('time','node');
nc{'wdir'}.long_name = 'Wave  Direction'; 
nc{'wdir'}.units     = 'degree'; 

nc{'tpeak'} = ncfloat('time','node');
nc{'tpeak'}.long_name = 'Relative Peak Period';
nc{'tpeak'}.units     = 's'; 

nc{'U10'} = ncfloat('time','node');
nc{'U10'}.long_name = 'Wind Velocity x-direction';
nc{'U10'}.units     = 'm/s'; 

nc{'V10'} = ncfloat('time','node');
nc{'V10'}.long_name = 'Wind Velocity y-direction';
nc{'V10'}.units     = 'm/s'; 

% static vars
nc{'x'}(:) = Xp;
nc{'y'}(:) = Yp;
nc{'h'}(:) = h;
nc{'nv'}(:,:) = tri';

% dump dynamic vars
for i=1:ntimes;

  fprintf('processing time %s\n',datestr(times(i)));
  %time
  shift = 678942.;  % datenum(2010,1,1,0,0,0)-greg2mjulian(2010,1,1,0,0,0);
  nc{'time'}(i) = times(i) - shift;

  %hs
  date = datestr(times(i),30);
  vname = ['Hsig_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  nc{'hs'}(i,:) = eval(vname)';

  %tp
  date = datestr(times(i),30);
  vname = ['RTpeak_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  nc{'tpeak'}(i,:) = eval(vname)';

  %depth
  date = datestr(times(i),30);
  vname = ['Depth_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  nc{'d'}(i,:) = eval(vname)';

  % wave dir
  date = datestr(times(i),30);
  vname = ['Dir_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  nc{'wdir'}(i,:) = eval(vname)';

  % U10 
  date = datestr(times(i),30);
  vname = ['Windv_x_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  var = eval(vname)';
  var(isnan(var)) = 0.0;
  nc{'U10'}(i,:) = var; 

  % V10 
  date = datestr(times(i),30);
  vname = ['Windv_y_',date(1:8),'_',date(10:15)];
  if(exist(vname) == 0)
    error('variable frame %s\n does not exist',vname)
  end;
  var = eval(vname)';
  var(isnan(var)) = 0.0;
  nc{'V10'}(i,:) = var'; 


end;


nc = close(nc);