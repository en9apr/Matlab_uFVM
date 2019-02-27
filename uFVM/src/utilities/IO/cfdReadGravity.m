function cfdReadGravity
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads g from "constant/g" file
%--------------------------------------------------------------------------
global Domain;

fprintf('\nReading Gravitational Properties ...\n');

projectDirectory = pwd;
projectDirectory = strrep(projectDirectory, '\', '/');

gFile = [projectDirectory, '/constant/g'];

% Check if "g" exists
if(~exist(gFile, 'file'))
    warning('\n%s\n', 'Assuming no gravitational effects');
    constant.name = 'g';
    constant.value = [0; 0; 0];
    cfdSetConstant(constant);
    return;
end

% Collect gravity attributes
[entry, dimensions] = getKeyValue('dimensions', gFile);
[entry, value] = getKeyValue('value', gFile);
%
% Store in foam data base
Domain.foam.g = {};
Domain.foam.g.dimensions = eval(dimensions{1});
Domain.foam.g.value = value{1};
%
% Create gravity mesh field
value = strtrim(Domain.foam.g.value);
C = textscan(value(2:end-1), '%f %f %f');
cfdSetupConstant('g', [C{1}; C{2}; C{3}]);
