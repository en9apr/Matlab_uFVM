function cfdReadSystem
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function read fvSchemes, fvSolution and controlDict files
%--------------------------------------------------------------------------
%
global Domain;

fprintf('\nReading System Dictionaries ...\n');
projectDirectory = pwd;
projectDirectory = strrep(projectDirectory, '\', '/');
 
% Read fvSchemes
ddtSchemes = readBlockData('ddtSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(ddtSchemes, 1)
    Domain.foam.fvSchemes.ddtSchemes{iScheme, 1} = ddtSchemes{iScheme, 1};
    Domain.foam.fvSchemes.ddtSchemes{iScheme, 2} = ddtSchemes{iScheme, 2};
end

gradSchemes = readBlockData('gradSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(gradSchemes, 1)
    Domain.foam.fvSchemes.gradSchemes{iScheme, 1} = gradSchemes{iScheme, 1};
    Domain.foam.fvSchemes.gradSchemes{iScheme, 2} = gradSchemes{iScheme, 2};
end

divSchemes = readBlockData('divSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(divSchemes, 1)
    Domain.foam.fvSchemes.divSchemes{iScheme, 1} = divSchemes{iScheme, 1};
    Domain.foam.fvSchemes.divSchemes{iScheme, 2} = divSchemes{iScheme, 2};
end

laplacianSchemes = readBlockData('laplacianSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(laplacianSchemes, 1)
    Domain.foam.fvSchemes.laplacianSchemes{iScheme, 1} = laplacianSchemes{iScheme, 1};
    Domain.foam.fvSchemes.laplacianSchemes{iScheme, 2} = laplacianSchemes{iScheme, 2};
end

interpolationSchemes = readBlockData('interpolationSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(interpolationSchemes, 1)
    Domain.foam.fvSchemes.interpolationSchemes{iScheme, 1} = interpolationSchemes{iScheme, 1};
    Domain.foam.fvSchemes.interpolationSchemes{iScheme, 2} = interpolationSchemes{iScheme, 2};
end

snGradSchemes = readBlockData('snGradSchemes', [projectDirectory, '/system/fvSchemes']);
for iScheme=1:size(snGradSchemes, 1)
    Domain.foam.fvSchemes.snGradSchemes{iScheme, 1} = snGradSchemes{iScheme, 1};
    Domain.foam.fvSchemes.snGradSchemes{iScheme, 2} = snGradSchemes{iScheme, 2};
end

% fluxRequired = readBlockData('fluxRequired', [projectDirectory, '/system/fvSchemes']);
% for iScheme=1:size(fluxRequired, 1)
%     Domain.foam.fvSchemes.fluxRequired{iScheme, 1} = fluxRequired{iScheme, 1};
%     Domain.foam.fvSchemes.fluxRequired{iScheme, 2} = fluxRequired{iScheme, 2};
% end

% Read controlDict
[entry, application] = getKeyValue('application', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{1, 1} = 'application';
Domain.foam.controlDict{1, 2} = application{1};

[entry, startFrom] = getKeyValue('startFrom', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{2, 1} = 'startFrom';
Domain.foam.controlDict{2, 2} = startFrom{1};

[entry, startTime] = getKeyValue('startTime', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{3, 1} = 'startTime';
Domain.foam.controlDict{3, 2} = str2double(startTime{1});

[entry, stopAt] = getKeyValue('stopAt', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{4, 1} = 'stopAt';
Domain.foam.controlDict{4, 2} = stopAt{1};

[entry, endTime] = getKeyValue('endTime', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{5, 1} = 'endTime';
Domain.foam.controlDict{5, 2} = str2double(endTime{1});

[entry, deltaT] = getKeyValue('deltaT', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{6, 1} = 'deltaT';
Domain.foam.controlDict{6, 2} = str2double(deltaT{1});

[entry, writeControl] = getKeyValue('writeControl', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{7, 1} = 'writeControl';
Domain.foam.controlDict{7, 2} = writeControl{1};

[entry, writeInterval] = getKeyValue('writeInterval', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{8, 1} = 'writeInterval';
Domain.foam.controlDict{8, 2} = str2double(writeInterval{1});

[entry, purgeWrite] = getKeyValue('purgeWrite', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{9, 1} = 'purgeWrite';
Domain.foam.controlDict{9, 2} = str2double(purgeWrite{1});

[entry, writeFormat] = getKeyValue('writeFormat', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{10, 1} = 'writeFormat';
Domain.foam.controlDict{10, 2} = writeFormat{1};

[entry, writePrecision] = getKeyValue('writePrecision', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{11, 1} = 'writePrecision';
Domain.foam.controlDict{11, 2} = str2double(writePrecision{1});

[entry, timeFormat] = getKeyValue('timeFormat', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{12, 1} = 'timeFormat';
Domain.foam.controlDict{12, 2} = timeFormat{1};

[entry, timePrecision] = getKeyValue('timePrecision', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{13, 1} = 'timePrecision';
Domain.foam.controlDict{13, 2} = str2double(timePrecision{1});

[entry, runTimeModifiable] = getKeyValue('runTimeModifiable', [projectDirectory, '/system/controlDict']);
Domain.foam.controlDict{14, 1} = 'runTimeModifiable';
Domain.foam.controlDict{14, 2} = runTimeModifiable{1};
%
% Read fvSolution
%
Domain.foam.fvSolution = {};

% Find algorithm
blockNames = getBlockNames([projectDirectory, '/system/fvSolution']);
for iBlockName=1:length(blockNames)
    if strcmp(blockNames{iBlockName}, 'SIMPLE')
        Domain.algorithm = 'SIMPLE';
        break;
    elseif strcmp(blockNames{iBlockName}, 'SIMPLEC')
        Domain.algorithm = 'SIMPLEC';
        break;
    elseif strcmp(blockNames{iBlockName}, 'SIMPLER')
        Domain.algorithm = 'SIMPLER';
        break;
    elseif strcmp(blockNames{iBlockName}, 'SIMPLEST')
        Domain.algorithm = 'SIMPLEST';
        break;
    elseif strcmp(blockNames{iBlockName}, 'PISO')
        Domain.algorithm = 'PISO';
        break;
    end
end

% Algorithm attributes
algorithm = Domain.algorithm;
algorithmData = readBlockData(algorithm, [projectDirectory, '/system/fvSolution']);

algorithmAttributes = {};

Domain.foam.fvSolution.(algorithm).nCorrectors = 1;
Domain.foam.fvSolution.(algorithm).nNonOrthogonalCorrectors = 0;

for iEntry=1:size(algorithmData,1)
    iMax = strfind(algorithmData{iEntry, 1}, 'Max');
    if ~isempty(iMax)
        fieldMax = getKeyValueFromBlock(algorithmData{iEntry, 1}, algorithm, [projectDirectory, '/system/fvSolution']);
        algorithmAttributes = setfield(algorithmAttributes, algorithmData{iEntry, 1}, fieldMax);
    end
    
    iMin = strfind(algorithmData{iEntry, 1}, 'Min');
    if ~isempty(iMin)
        fieldMin = getKeyValueFromBlock(algorithmData{iEntry, 1}, algorithm, [projectDirectory, '/system/fvSolution']);
        algorithmAttributes = setfield(algorithmAttributes, algorithmData{iEntry, 1}, fieldMin);
    end 
        
    if strcmp(algorithmData{iEntry, 1}, 'nCorrectors') || strcmp(algorithmData{iEntry, 1}, 'nNonOrthogonalCorrectors')
        if ~isempty(str2double(algorithmData{iEntry, 2}))
            algorithmAttributes = setfield(algorithmAttributes, algorithmData{iEntry, 1}, str2double(algorithmData{iEntry, 2}));
        end
    end 
    
    % Reference cell and value for pressure
    if strcmp(algorithmData{iEntry, 1}, 'pRefCell')
        if ~isempty(str2double(algorithmData{iEntry, 2}))
            if str2double(algorithmData{iEntry, 2})==0
                index = 1;
            else
                index = 0;
            end
            algorithmAttributes = setfield(algorithmAttributes, algorithmData{iEntry, 1}, str2double(algorithmData{iEntry, 2})+index);
        end
    end 
    if strcmp(algorithmData{iEntry, 1}, 'pRefValue')
        if ~isempty(str2double(algorithmData{iEntry, 2}))
            algorithmAttributes = setfield(algorithmAttributes, algorithmData{iEntry, 1}, str2double(algorithmData{iEntry, 2}));
        end
    end     
end

Domain.foam.fvSolution = setfield(Domain.foam.fvSolution, algorithm, algorithmAttributes);

% Set default values for nCorrectors and nNonOrthogonalCorrectors
if ~isfield(Domain.foam.fvSolution.(algorithm), 'nCorrectors')
    Domain.foam.fvSolution.(algorithm).nCorrectors = 1;
end
if ~isfield(Domain.foam.fvSolution.(algorithm), 'nNonOrthogonalCorrectors')
    Domain.foam.fvSolution.(algorithm).nNonOrthogonalCorrectors = 0;
end

% Linear Solver Settings
% Read files from 0 directory
files = dir([projectDirectory, '/0/']);

Domain.foam.fvSolution.solvers = {};

iField = 1;
for file=files'
    fieldName = file.name;    
    
    if file.bytes==0 || file.isdir        
        continue;
    end
    
    fieldSolution = readBlock(['solvers/',fieldName], [projectDirectory, '/system/fvSolution']);
    
    if ~isempty(fieldSolution)
        Domain.foam.fvSolution.solvers = setfield(Domain.foam.fvSolution.solvers, cfdConvertFieldPhaseName(fieldName), fieldSolution);
        
        % Residual Tolerance
        Domain.foam.fvSolution.(algorithm).residualControl{iField, 1} = cfdConvertFieldPhaseName(fieldName);
        fieldResidualControl = getKeyValueFromBlock(fieldName, 'residualControl', [projectDirectory, '/system/fvSolution']);
        if isempty(fieldResidualControl)
            Domain.foam.fvSolution.(algorithm).residualControl{iField, 2} = '1e-6';
        else
            Domain.foam.fvSolution.(algorithm).residualControl{iField, 2} = getKeyValueFromBlock(fieldName, 'residualControl', [projectDirectory, '/system/fvSolution']);
        end
        iField = iField + 1;
    end
end

% Relaxation Factors
fields = readBlock('relaxationFactors/fields', [projectDirectory, '/system/fvSolution']);
Domain.foam.fvSolution.relaxationFactors.fields = fields;

equations = readBlock('relaxationFactors/equations', [projectDirectory, '/system/fvSolution']);
Domain.foam.fvSolution.relaxationFactors.equations = equations;

% Additional:
% Check if all physical types of bc's are wall
cavity = true;
theMesh = cfdGetMesh;
for iBoundary=1:length(theMesh.boundaries)
    boundaryType = theMesh.boundaries(iBoundary).type;
    if ~strcmp(boundaryType, 'wall') && ~strcmp(boundaryType, 'empty') && ~strcmp(boundaryType, 'symmetry')
        cavity = false;
        break;
    end
end

applicationClass = cfdGetApplicationClass;

if cavity && ~strcmp(applicationClass, 'basic')
   cfdSetFixedElement;
end



