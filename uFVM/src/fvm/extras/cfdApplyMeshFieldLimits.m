function theMeshField = cfdApplyMeshFieldLimits(theMeshField)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function crops the mesh field values at any element if the
%   corresponding value lies outside the limits. Currently working for
%   scalar fields
%--------------------------------------------------------------------------

global Domain;
algorithm = cfdGetAlgorithm;

if isempty(Domain.foam.fvSolution)
    return;
end

algorithmAttributes = Domain.foam.fvSolution.(algorithm);
attributeNames = fieldnames(algorithmAttributes);

if any(ismember(attributeNames, [theMeshField.name, 'Max']))
    if strcmp(theMeshField.type, 'Scalar')
        theMeshField.phi = min(theMeshField.phi, str2double(algorithmAttributes.([theMeshField.name, 'Max'])));        
    end        
end

if any(ismember(attributeNames, [theMeshField.name, 'Min']))
    if strcmp(theMeshField.type, 'Scalar')
        theMeshField.phi = max(theMeshField.phi, str2double(algorithmAttributes.([theMeshField.name, 'Min'])));
    end        
end
