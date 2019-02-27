function cfdDefineApplicationClass(applicationClass)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets the class of application. The class could be
%   incompressible, compressible, multiphase, heatTransfer
%--------------------------------------------------------------------------
registeredApplicationClasses = {'basic', 'incompressible', 'compressible', 'multiphase', 'heatTransfer'};
if(~isempty(strcmp(applicationClass, registeredApplicationClasses)))
    cfdSetApplicationClass(applicationClass);
else
    error('The application class entered can''t be recognized, choose one of the following:\n%s\n%s\n%s\n%s', ...
        'incompressible', 'compressible', 'multiphase', 'heatTransfer');    
end