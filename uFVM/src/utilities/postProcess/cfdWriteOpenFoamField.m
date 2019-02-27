function cfdWriteOpenFoamField(theField, timeStep)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function writes the field in OpenFOAM format
%--------------------------------------------------------------------------

global Domain;
theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBElements = theMesh.numberOfBElements;

theFieldFoamAttributes = Domain.foam.fields.(theField.name);
theFieldName = cfdConvertFieldPhaseName(theField.name);

% Concatenate header
text = ...
    {'/*--------------------------------*- C++ -*----------------------------------*\\', ...
    '| =========                 |                                                 |', ...
    '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |', ...
    '|  \\    /   O peration     | Version:  v1706                                 |', ...
    '|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |', ...
    '|    \\/     M anipulation  |                                                 |', ...
    '\\*---------------------------------------------------------------------------*/', ...
    'FoamFile', ...
    '{', ...
    '    version     2.0;', ...
    '    format      ascii;', ...
    ['    class       ', theFieldFoamAttributes.class,';'], ...
    ['    location    "', num2str(timeStep),'";'], ...
    ['    object      ', theFieldName,';'], ...
    '}', ...
    '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //', ...
    ''};

% Concatenate dimensions
text{end+1} = ['dimensions      ', theFieldFoamAttributes.dimensions];
text{end+1} = '';
text{end+1} = '';

% Concatenate internal field
if strcmp(theFieldFoamAttributes.class, 'volScalarField')
    text{end+1} = 'internalField   nonuniform List<scalar>';
    text{end+1} = num2str(numberOfElements);
    text{end+1} = '(';
    for iElement=1:numberOfElements
        text{end+1} = num2str(theField.phi(iElement));
    end
    text{end+1} = ')';
    text{end+1} = ';';
elseif strcmp(theFieldFoamAttributes.class, 'volVectorField')
    text{end+1} = 'internalField   nonuniform List<vector>';
    text{end+1} = num2str(numberOfElements);
    text{end+1} = '(';
    for iElement=1:numberOfElements
        text{end+1} = ['(', strjoin({num2str(theField.phi(iElement, 1)), num2str(theField.phi(iElement, 2)), num2str(theField.phi(iElement, 3))}), ')'];
    end
    text{end+1} = ')';
    text{end+1} = ';';
elseif strcmp(theFieldFoamAttributes.class, 'surfaceScalarField')
    text{end+1} = 'internalField   nonuniform List<scalar>';
    text{end+1} = num2str(numberOfInteriorFaces);
    text{end+1} = '(';
    for iFace=1:numberOfInteriorFaces
        text{end+1} = num2str(theField.phi(iFace));
    end
    text{end+1} = ')';
    text{end+1} = ';';
end


% Concatenate boundary field
text{end+1} = '';
text{end+1} = 'boundaryField';
text{end+1} = '{';

theNumberOfPatches = theMesh.numberOfBoundaries;
for iPatch=1:theNumberOfPatches
    theBoundary = theMesh.boundaries(iPatch);
    numberOfBFaces = theBoundary.numberOfBFaces;
    
    %
    iFaceStart = theBoundary.startFace;
    iFaceEnd = iFaceStart+numberOfBFaces-1;
    iBFaces = iFaceStart:iFaceEnd;
    %
    iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
    iElementEnd = iElementStart+numberOfBFaces-1;
    iBElements = iElementStart:iElementEnd;
    
    text{end+1} = ['    ', theBoundary.userName];
    text{end+1} = '    {';
    
    boundaryField = theFieldFoamAttributes.boundaryField{iPatch};
    type = boundaryField.type;
    text{end+1} = ['        type             ', type, ';'];
    
    if isfield(boundaryField, 'value')
        if strcmp(boundaryField.valueType, 'nonuniform')
            if strcmp(theFieldFoamAttributes.class, 'volScalarField')
                text{end+1} = '        value            nonuniform List<scalar>';
                text{end+1} = num2str(numberOfBFaces);
                text{end+1} = '(';
                for iBElement=iElementStart:iElementEnd
                    text{end+1} = num2str(theField.phi(iBElement));
                end
                text{end+1} = ')';
                text{end+1} = ';';
            elseif strcmp(theFieldFoamAttributes.class, 'surfaceScalarField')
                text{end+1} = '        value            nonuniform List<scalar>';
                text{end+1} = num2str(numberOfBFaces);
                text{end+1} = '(';
                for iBFace=iFaceStart:iFaceEnd
                    text{end+1} = num2str(theField.phi(iBFace));
                end
                text{end+1} = ')';
                text{end+1} = ';';                
            elseif strcmp(theFieldFoamAttributes.class, 'volVectorField')
                text{end+1} = '        value            nonuniform List<vector>';
                text{end+1} = num2str(numberOfBFaces);
                text{end+1} = '(';
                for iBElement=iElementStart:iElementEnd
                    text{end+1} = ['(', strjoin({num2str(theField.phi(iBElement, 1)), num2str(theField.phi(iBElement, 2)), num2str(theField.phi(iBElement, 3))}), ')'];
                end
                text{end+1} = ')';
                text{end+1} = ';';
            end
        else
            if strcmp(theFieldFoamAttributes.class, 'volScalarField') || strcmp(theFieldFoamAttributes.class, 'surfaceScalarField')
                text{end+1} = ['        value            uniform ', num2str(boundaryField.value), ';'];
            else
                text{end+1} = ['        value            uniform (', strjoin({num2str(boundaryField.value(1)), num2str(boundaryField.value(2)), num2str(boundaryField.value(3))}), ');'];
            end            
        end
    end
    
    text{end+1} = '    }';
    
end
text{end+1} = '}';
text{end+1} = '';
text{end+1} = '// ************************************************************************* //';

% Print to time directory
fileID = fopen([num2str(timeStep), '\', theFieldName], 'w');
for i=1:length(text)
    fprintf(fileID, '%s\n', text{i});
end
fclose(fileID);
