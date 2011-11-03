function set_spectide(Mobj)

% Setup spectral tides on the open boundary and dump a spectral file  
%
% function set_spectide(Mobj)  
%
% DESCRIPTION:
%    Setup spectral tides on the open boundary and dump a spectral file
%    This is a USER DEFINED driver program for the FVCOM spectral tide
%    It requires USER Modification to work 
%
% INPUT
%    Mobj         = Matlab mesh object
%
% OUTPUT:
%
% EXAMPLE USAGE
%    set_spectide(Mobj)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'set_spectide';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Set Constants
%------------------------------------------------------------------------------
MyTitle = 'spectral tide test';
SpectralFile = 'tst_spectide.nc';

%------------------------------------------------------------------------------
% Set Component Periods
%------------------------------------------------------------------------------
Components = {   'M2',    'N2',    'S2',   'O1',    'K1'};
Period     = [44714.16, 45570.05, 43200, 92949.63, 86164.09];

%------------------------------------------------------------------------------
% Setup user defined phase and amplitude along the open boundaries
% need to set:
%   1.) Period - vector containing component period in seconds
%   2.) Amp    - array of size [Nobcs, Ncomponents] containing amplitude in m
%   3.) Phase  - array of size [Nobcs, Ncomponents] containing phase in degrees
%------------------------------------------------------------------------------
nComps = 1;

if(Mobj.nObs==0)
	warning('cannot setup spectral open boundary, there is no open boundary in the mesh struct')
	return
end;

cnt = 0;
for ob=1:Mobj.nObs
	nObcs = Mobj.nObcNodes(ob);
	for j=1:nObcs
		cnt = cnt + 1;
		ObcNodes(cnt) = Mobj.obc_nodes(1,j);  %set the open boundary nodes
		for i=1:nObcs
			if(ob==1)
				Amp(cnt,1:nComps) = [1.0];  
			else
				Amp(cnt,1:nComps) = [1.0];
			end;
		end;

		for i=1:nObcs
			Phase(cnt,1:nComps) = [0.];
		end;
	end;
end;

%------------------------------------------------------------------------------
% Dump a spectral tide file in NetCDF
%------------------------------------------------------------------------------
write_FVCOM_spectide(ObcNodes,Period(1:nComps),Phase,Amp,SpectralFile,MyTitle)

if(ftbverbose); fprintf(['end   : ' subname '\n']);end;
