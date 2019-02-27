function cfdReadTurbulenceProperties
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads turbulence properties from
%   "constant/turbulenceProperties" file
%--------------------------------------------------------------------------

global Domain;

fprintf('\n\nReading Turbulence Properties ...\n');

projectDirectory = pwd;
projectDirectory = strrep(projectDirectory, '\', '/');

turbulencePropertiesFile = [projectDirectory, '/constant/turbulenceProperties'];

% Check if "transportProperties" exists
if(~exist(turbulencePropertiesFile, 'file'))
    return;
%     error('\n%s\n','"turbulenceProperties" file is not found in "~foamDirectory/constant"');
end

% get simulation type
[entry, simulationType] = getKeyValue('simulationType', turbulencePropertiesFile);
simulationTypeDetails = readBlockData(simulationType, turbulencePropertiesFile);
for iDetail=1:size(simulationTypeDetails, 1)
    Domain.foam.turbulenceProperties{iDetail, 1} = simulationTypeDetails{iDetail, 1};
    Domain.foam.turbulenceProperties{iDetail, 2} = simulationTypeDetails{iDetail, 2};
end
