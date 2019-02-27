function cfdReadOpenFoamBCs(caseDirectory)
%===================================================

%  written by the CFD Group @ AUB, Fall 2017
%===================================================
if(nargin==0)
    caseDirectory = '~';
    caseDirectory = uigetdir(caseDirectory);
end

fprintf('\nReading OpenFOAM System Dictionaries ...\n\n');

fvSchemesFile = [caseDirectory,'/system/fvSchemes'];
controlDictFile = [caseDirectory,'/system/controlDict'];
fvSolutionFile = [caseDirectory,'/system/fvSolution'];

%
% Read time settings
%
ddtSchemes_settings = cfdReadOpenFoamBlockFromFile(fvSchemesFile, 'ddtSchemes');
startTime = cfdReadOpenFoamKeyValue(controlDictFile, 'startTime');
endTime = cfdReadOpenFoamKeyValue(controlDictFile, 'endTime');
dt = cfdReadOpenFoamKeyValue(controlDictFile, 'deltaT');

if strcmp(ddtSchemes_settings{2}, 'steadyState')
    cfdSetupTime('type', 'Steady', 'startTime', startTime, 'endTime', endTime, 'dt', dt);
elseif strcmp(ddtSchemes_settings{2}, 'Euler')
    cfdSetupTime('type', 'Transient', 'startTime', startTime, 'endTime', endTime, 'dt', dt);
end