function m_mappath(TT)
% Sets up path for codas software
% one optional parameter is allowed, empty for adding the path
% 0 or anyother thing for removing it.
 % addpath d:/matlab/codas
if (~nargin) 
   addpath /users/modellers/rito/matlab/m_map
else
   rmpath /users/modellers/rito/matlab/m_map
end

   