function theInitialField = cfdGenerateInitialResults(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads initial equation field and creates the
%   corresponding mesh field
%--------------------------------------------------------------------------

global Domain;
projectDirectory = cd;
if isempty(Domain) || ~isfield(Domain, 'projectDirectory')
    % Check if mesh info available already in constant/polyMesh
    meshExists = exist('constant/polyMesh/theMesh.mat', 'file');
    if meshExists==2
        load([projectDirectory, '\constant\polyMesh\theMesh.mat']);
        cfdSetMesh(theMesh);
    else
        cfdReadPolyMesh;
        theMesh = cfdGetMesh;
        save([projectDirectory, '\constant\polyMesh\theMesh.mat'], 'theMesh');
    end    
    
    % Read fields from OpenFOAM files
    cfdReadFields;
else
    theMesh = cfdGetMesh;
end

% Get field information from data base
theField = getfield(Domain.foam.fields, theEquationName);

% The Equation Initial Conditions
theFieldIC = theField.internalField;
if strfind(theFieldIC, 'uniform')
    if strcmp(theField.class, 'volVectorField')
        op = strfind(theField.internalFieldValue, '(');
        cp = strfind(theField.internalFieldValue, ')');
        C = textscan(theField.internalFieldValue(op+1:cp-1), '%f %f %f');
        ic = ['[',num2str(C{1}),';',num2str(C{2}), ...
            ';',num2str(C{3}),']'];
    else
        C = textscan(theField.internalFieldValue, '%s %f');
        ic = num2str(C{2});
    end
else
    ic = theField.internalFieldValue;
end

if strcmp(theField.class,'volScalarField')
    theInitialField.type = 'Scalar';
elseif strcmp(theField.class,'volVectorField')
    theInitialField.type = 'Vector';
end

theNumberOfElements = theMesh.numberOfElements;
theNumberOfInteriorFaces = theMesh.numberOfInteriorFaces;
theNumberOfFaces = theMesh.numberOfFaces;

theNumberOfBElements =  theMesh.numberOfBElements;
theElementCentroids = [theMesh.elements(:).centroid]';
theBElementCentroids = [theMesh.faces(theNumberOfInteriorFaces+1:theNumberOfFaces).centroid]';

theCentroids = [theElementCentroids; theBElementCentroids];

x = theCentroids(:,1);
y = theCentroids(:,2);
z = theCentroids(:,3);

if strcmp(theInitialField.type,'Scalar')
    theInitialField.phi = eval(ic) .* ones(theNumberOfElements+theNumberOfBElements,1);
elseif strcmp(theInitialField.type,'Vector')
    for iComponent=1:3
        theComponentFormula = cfdGetFormulaForComponent(ic,iComponent);
        theInitialField.phi(:,iComponent) = eval(theComponentFormula) .* ones(theNumberOfElements+theNumberOfBElements,1);
    end
end



