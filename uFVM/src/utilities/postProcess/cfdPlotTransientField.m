function cfdPlotTransientField(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function plots the contour of the field in transient mode
%--------------------------------------------------------------------------

% Re-construct domain if not available
cfdReconstructDomain;

figure;
axis off
axis equal

% Get mesh
theMesh = cfdGetMesh;

theNumberOfNodes = theMesh.numberOfNodes;
theNumberOfFaces = theMesh.numberOfFaces;

for iNode=1:theNumberOfNodes
    v(iNode,:) = theMesh.nodes(iNode).centroid';
end

for iFace=1:theNumberOfFaces
    f{iFace,:} = theMesh.faces(iFace).iNodes;
end

f = cfdConvertCelltoMatrix(f);

% Get min and max phi
timeSteps = cfdGetTimeSteps;

phiMax = zeros(length(timeSteps), 1);
phiMin = zeros(length(timeSteps), 1);
i = 1;
for timeStep=timeSteps'
    theFoamFields = readTimeDirectory(num2str(timeStep), theFieldName);
    theFoamField = getfield(theFoamFields, theFieldName);
    theMeshField = createMeshField(theFieldName, theFoamField);
    
    if strcmp(theMeshField.type, 'Scalar')
        phiMax(i, 1) = max(theMeshField.phi);
        phiMin(i, 1) = min(theMeshField.phi);
    else
        phiMax(i, 1) = max(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
        phiMax(i, 1) = min(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
    end
    i = i + 1;
end

phiMax = max(phiMax);
phiMin = min(phiMin);

% Plot
for timeStep=timeSteps'
    if timeStep==0
        theFoamFields = readTimeDirectory(num2str(timeStep), theFieldName);
        theFoamField = getfield(theFoamFields, theFieldName);
        theMeshField = createMeshField(theFieldName, theFoamField);
        
        if strcmp(theMeshField.type, 'Vector')
            phiNodes = cfdInterpolateFromElementsToNodes(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
        else
            phiNodes = cfdInterpolateFromElementsToNodes(theMeshField.phi(:,1));
        end
        p = patch('vertices', v, 'faces', f, 'FaceVertexCData', phiNodes, 'FaceColor', 'interp', 'EdgeAlpha', .1);
        
        caxis([phiMin, phiMax]);
        
        rotate3d;
        set(gca,'Clipping','off');
        view(30,40)
        
        colorbar;
        colormap jet;
        title(theFieldName);
    else
        theFoamFields = readTimeDirectory(num2str(timeStep), theFieldName);
        theFoamField = getfield(theFoamFields, theFieldName);
        theMeshField = createMeshField(theFieldName, theFoamField);
        
        if strcmp(theMeshField.type, 'Vector')
            phiNodes = cfdInterpolateFromElementsToNodes(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
        else
            phiNodes = cfdInterpolateFromElementsToNodes(theMeshField.phi(:,1));
        end
        set(p, 'CData', phiNodes);
        drawnow;
        
    end
end