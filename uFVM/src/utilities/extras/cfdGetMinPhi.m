function phiMin = cfdGetMinPhi(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function collects the minimum phi from all results
%--------------------------------------------------------------------------

timeSteps = cfdGetTimeSteps;
phiMin = zeros(length(timeSteps), 1);
i = 1;
for timeStep=timeSteps'
    theFoamFields = readTimeDirectory(num2str(timeStep));
    theFoamField = getfield(theFoamFields, theFieldName);
    theMeshField = createMeshField(theFieldName, theFoamField);
    
    if strcmp(theMeshField.type, 'Scalar')
        phiMin(i, 1) = min(theMeshField.phi);
    else
        phiMin(i, 1) = min(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
    end
    i = i + 1;
end

phiMin = min(phiMin);