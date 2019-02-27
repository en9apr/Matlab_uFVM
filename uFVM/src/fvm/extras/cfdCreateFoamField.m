function field = cfdCreateFoamField(object, class, dimensions)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function creates a field in the OpenFOAM format
%--------------------------------------------------------------------------

% Create corresponding field in Foam format
theMesh = cfdGetMesh;

% Create the foam field which is to be used
field = {};
field.class = class;
field.dimensions = dimensions;

if strcmp(object, 'mdot_f')
    % Create the internal field
    
    % Get mesh details
    numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
    
    field.internalField.valueType = 'nonuniform';
    field.internalField.value = zeros(numberOfInteriorFaces, 1);
    
    % Create the boundary field
    boundaries = theMesh.boundaries;
    
    for iBoundary=1:length(boundaries)
        % Fix calculated boundary type for all patches
        field.boundaryField{iBoundary}.type = 'calculated';
        
        numberOfBFaces = boundaries(iBoundary).numberOfBFaces;
        
        thePhysicalType = boundaries(iBoundary).type;
        if strcmp(thePhysicalType,'wall') || strcmp(thePhysicalType,'symmetry') || strcmp(thePhysicalType,'empty')
            field.boundaryField{iBoundary}.valueType = 'uniform';
            field.boundaryField{iBoundary}.value = 0;
        elseif strcmp(thePhysicalType,'inlet') || strcmp(thePhysicalType,'outlet')
            field.boundaryField{iBoundary}.valueType = 'nonuniform';
            field.boundaryField{iBoundary}.value = zeros(numberOfBFaces, 1);
        end
    end        
else
    
    
    
end