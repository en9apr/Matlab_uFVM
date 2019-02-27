function cfdReconstructDomain
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function re-constructs the domain by reading the existing data as
%   well latest results
%--------------------------------------------------------------------------

% Check if Domain exists
global Domain;

if ~isempty(Domain)
    return;
end

% Create and store mesh
theMesh = readPolyMesh;
cfdSetMesh(theMesh);

% Read fields
timeSteps = cfdGetTimeSteps;
latestTimeStep = max(timeSteps);

theFoamFields = readTimeDirectory(num2str(latestTimeStep));
theFieldNames = fieldnames(theFoamFields);

Domain.Step0.Elements.fields = {};
Domain.Step0.Faces.fields = {};
Domain.Step0.Nodes.fields = {};

for iField=1:length(theFieldNames)
    theFoamField = getfield(theFoamFields, theFieldNames{iField});
    theMeshField = createMeshField(theFieldNames{iField}, theFoamField);
    cfdSetMeshField(theMeshField);    
end

