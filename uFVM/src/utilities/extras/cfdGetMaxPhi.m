function phiMax = cfdGetMaxPhi(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function collects the maximum phi from all results
%--------------------------------------------------------------------------

timeSteps = cfdGetTimeSteps;
phiMax = zeros(length(timeSteps), 1);
i = 1;
for timeStep=timeSteps'
    theFoamFields = readTimeDirectory(num2str(timeStep));
    theFoamField = getfield(theFoamFields, theFieldName);
    theMeshField = createMeshField(theFieldName, theFoamField);
    
    if strcmp(theMeshField.type, 'Scalar')
        phiMax(i, 1) = max(theMeshField.phi);
    else
        phiMax(i, 1) = max(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
    end
    i = i + 1;
end

phiMax = max(phiMax);