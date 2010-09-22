function write_FVCOM_spectide(ObcNodes,Period,Phase,Amp,SpectralFile,MyTitle)   
	
% Write an FVCOM spectral tidal elevation forcing file 
%
% function write_FVCOM_spectide(ObcNodes,Period,Phase,Amp,SpectralFile,MyTitle)   
%
% DESCRIPTION:
%    Write an FVCOM NetCDF spectral tidal elevation forcing file
%
% INPUT:
%   ObcNodes     = list of open boundary nodes of size [nObcs]
%   Period       = list of periods in seconds of size [nComponents]
%   Phase        = list of phases in degrees of size [nObcs,nComponents]
%   Amp          = list of amplitudes (m) of size [nObcs,nComponents]
%   SpectralFile = name of NetCDF file
%   MyTitle      = case title, written as global attribute of NetCDF file
%
% OUTPUT:
%    SpectralFile, A NetCDF FVCOM spectral tide forcing file
%
% EXAMPLE USAGE
%    write_FVCOM_spectide(ObcNodes,Period,Phase,Amp,SpectralFile,MyTitle) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
warning off

subname = 'write_FVCOM_spectide';
fprintf('\n')
fprintf(['begin : ' subname '\n'])
%------------------------------------------------------------------------------
% Sanity check on input and dimensions
%------------------------------------------------------------------------------
nComponents = prod(size(Period));
fprintf('Number of Tide Components %d\n',nComponents)

nObcs = prod(size(ObcNodes));
fprintf('Number of Open Boundary Nodes %d\n',nObcs)

[chk1,chk2] = size(Amp);
if( (nObcs-chk1)*(nComponents-chk2) ~= 0)
	fprintf('Amp dimensions do not match!!!')
	fprintf('nObcs %d nComponens %d\n',chk1,chk2)
	stop
end;

[chk1,chk2] = size(Phase);
if( (nObcs-chk1)*(nComponents-chk2) ~= 0)
	fprintf('Phase dimensions do not match!!!')
	fprintf('nObcs %d nComponens %d\n',chk1,chk2)
	stop
end;

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------

nc = netcdf(SpectralFile, 'clobber');             

% dump header
nc.type = 'FVCOM SPECTRAL ELEVATION FORCING FILE' ;
nc.title = MyTitle;
nc.history = 'FILE CREATED using write_FVCOM_spectide';

% dimensions
nc('nobc') = nObcs;
nc('tidal_components') = nComponents;
nc('DateStrLen') = 26;

% variables
nc{'obc_nodes'} = ncint('nobc');
nc{'obc_nodes'}.long_name = 'Open Boundary Node Number'; 
nc{'obc_nodes'}.grid = 'obc_grid'; 

nc{'tide_period'} = ncfloat('tidal_components');
nc{'tide_period'}.long_name = 'tide angular period';
nc{'tide_period'}.units = 'seconds';

nc{'tide_Eref'} = ncfloat('nobc');
nc{'tide_Eref'}.long_name = 'tidal elevation reference level';
nc{'tide_Eref'}.units = 'meters';

nc{'tide_Ephase'} = ncfloat('tidal_components', 'nobc');
nc{'tide_Ephase'}.long_name = 'tidal elevation phase angle';
nc{'tide_Ephase'}.units = 'degrees, time of maximum elevation with respect to chosen time origin';

nc{'tide_Eamp'} = ncfloat('tidal_components', 'nobc');
nc{'tide_Eamp'}.long_name = 'tidal elevation amplitude';
nc{'tide_Eamp'}.units = 'meters';

nc{'equilibrium_tide_Eamp'} = ncfloat('tidal_components');
nc{'equilibrium_tide_Eamp'}.long_name = 'equilibrium tidal elevation amplitude';
nc{'equilibrium_tide_Eamp'}.units = 'meters';

nc{'equilibrium_beta_love'} = ncfloat('tidal_components');
nc{'equilibrium_beta_love'}.formula = 'beta=1+klove-hlove';

nc{'equilibrium_tide_type'} = ncchar('tidal_components', 'DateStrLen');
nc{'equilibrium_tide_type'}.long_name = 'formula';
nc{'equilibrium_tide_type'}.units = 'beta=1+klove-hlove';

nc{'time_origin'} = ncfloat;
nc{'time_origin'}.long_name = 'time';
nc{'time_origin'}.units = 'days since 0.0';
nc{'time_origin'}.time_zone = 'none';


% data
nc{'obc_nodes'}(:)   = ObcNodes; 
nc{'tide_period'}(:) = Period;
nc{'tide_Eref'}(:)   = 0.; 
nc{'tide_Ephase'}(:,:) = Phase';
nc{'tide_Eamp'}(:,:)   = Amp';
nc{'equilibrium_tide_Eamp'}(:,:) = 0.0;
nc{'equilibrium_beta_love'}(:,:) = 0.0;

for i=1:nComponents
  nc{'equilibrium_tide_type'}(i,1:26) = 'SEMIDIURNAL               ';
end;
nc{'time_origin'}(:) = 0.0;

nc = close(nc);    


fprintf(['end   : ' subname '\n'])

