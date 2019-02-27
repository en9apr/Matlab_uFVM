function fields = readTimeDirectory(timeDirectory, theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads fields from time directory
%--------------------------------------------------------------------------

% Re-construct domain if not available
cfdReconstructDomain;

% Read field
files = dir(timeDirectory);

fields = {};
projectDirectory = pwd;
if strfind(projectDirectory, '/')
    projectDirectory = strrep(projectDirectory, '/', '\');
end

theMesh = cfdGetMesh;

for iFile=1:length(files)
    if (files(iFile).bytes)==0 || (files(iFile).isdir)
        continue;
    end
    
    % get field name from file name
    fieldName = files(iFile).name;
    
    if nargin==2
        if ~strcmp(fieldName, theFieldName)
            continue;
        end
    end
    
    % Initialize field structure to store the data inside
    field = {};
    
    % Read class (volScalarField, volVectorField)
    [entry, class] = getKeyValue('class', [projectDirectory, '\', timeDirectory,'\', fieldName], false);
    field.class = class{1};
    
    % Read dimensions
    [entry, dimensions] = getKeyValue('dimensions', [projectDirectory, '\', timeDirectory,'\', fieldName]);
    field.dimensions = dimensions{1};
    
    % Read internal field
    [entry, internalField] = getKeyValue('internalField', [projectDirectory, '\', timeDirectory,'\', fieldName]);
    C = textscan(internalField{1}, '%s', 1);
    if strcmp(C{1}{1}, 'uniform')
        field.internalField.valueType = 'uniform';
        if strcmp(field.class, 'volScalarField') || strcmp(field.class, 'surfaceScalarField')
            value_str = textscan(internalField{1}, 'uniform %f;', 1);
            field.internalField.value = value_str{1};
        elseif strcmp(field.class, 'volVectorField')
            value_str = textscan(internalField{1}, 'uniform (%f %f %f);', 1);
            field.internalField.value = [value_str{1}, value_str{2}, value_str{3}];
        end
    elseif strcmp(C{1}{1}, 'nonuniform')
        field.internalField.valueType = 'nonuniform';
        list = readNonuniformList('internalField', [projectDirectory, '\', timeDirectory,'\', fieldName]);
        field.internalField.value = list;
    end
    
    % Read boundary field
    boundaries = theMesh.boundaries;
    field.boundaryField = {};
    
    for iBoundary=1:length(boundaries)
        type = getKeyValueFromBlock('type', ['boundaryField/', boundaries(iBoundary).userName], [projectDirectory, '\', timeDirectory,'\', fieldName]);
        if strcmp(type, 'calculated') || strcmp(type, 'fixedValue') || strcmp(type, 'fixedGradient')
            field.boundaryField{iBoundary}.type = type;
            value = getKeyValueFromBlock('value', ['boundaryField/', boundaries(iBoundary).userName], [projectDirectory, '\', timeDirectory,'\', fieldName]);
            C = textscan(value, '%s', 1);
            if strcmp(C{1}{1}, 'uniform')
                field.boundaryField{iBoundary}.valueType = 'uniform';
                if strcmp(field.class, 'volScalarField') || strcmp(field.class, 'surfaceScalarField')
                    value_str = textscan(value, 'uniform %f', 1);
                    field.boundaryField{iBoundary}.value = value_str{1};
                elseif strcmp(field.class, 'volVectorField')
                    value_str = textscan(value, 'uniform (%f %f %f)', 1);
                    field.boundaryField{iBoundary}.value = [value_str{1}, value_str{2}, value_str{3}];
                end
            elseif strcmp(C{1}{1}, 'nonuniform')
                list = readNonuniformList(boundaries(iBoundary).userName, [projectDirectory, '\', timeDirectory,'\', fieldName]);
                field.boundaryField{iBoundary}.valueType = 'nonuniform';
                field.boundaryField{iBoundary}.value = list;
            end
        elseif strcmp(type, 'zeroGradient') || strcmp(type, 'noSlip') || strcmp(type, 'slip') || strcmp(type, 'empty') || strcmp(type, 'symmetry')
            field.boundaryField{iBoundary}.type = type;
        else
            boundaryField = readBlockData(boundaries(iBoundary).userName, [projectDirectory, '\', timeDirectory,'\', fieldName]);
            field.boundaryField{iBoundary} = {};
            for iBAttribute=1:size(boundaryField, 1)
                field.boundaryField{iBoundary} = setfield(field.boundaryField{iBoundary}, boundaryField{iBAttribute, 1}, boundaryField{iBAttribute, 2});
            end
        end
    end
    
    fields = setfield(fields, cfdConvertFieldPhaseName(fieldName), field);
end