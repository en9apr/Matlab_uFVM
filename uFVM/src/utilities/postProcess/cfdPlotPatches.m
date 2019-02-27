function cfdPlotPatches(iPatches)

figure;

theMesh = cfdGetMesh;

axis off
axis equal

theNumberOfNodes = theMesh.numberOfNodes;
theNumberOfFaces = theMesh.numberOfFaces;

for iNode=1:theNumberOfNodes
    v(iNode,:) = theMesh.nodes(iNode).centroid';
end

for iFace=1:theNumberOfFaces
    f{iFace,:} = theMesh.faces(iFace).iNodes;
end

f = cfdConvertCelltoMatrix(f);

phiNodes = zeros(length(v),1);
patch('vertices', v, 'faces', f, 'FaceVertexCData', phiNodes,'FaceAlpha', 0.1, 'EdgeAlpha', .4);

for iPatch=iPatches
    if((iPatch<=theMesh.numberOfBoundaries) && (iPatch>0))
        theBoundary = theMesh.boundaries(iPatch);
        iFaces = theBoundary.startFace:theBoundary.startFace+theBoundary.numberOfBFaces-1;
        for iFace=iFaces
            theFace = theMesh.faces(iFace);
            theNodes = theMesh.nodes(theFace.iNodes);
            local_phi = iPatch;
            XYZ = [theNodes.centroid]';
            patch(XYZ(:,1), XYZ(:,2), XYZ(:,3), local_phi, 'FaceAlpha', 1, 'EdgeAlpha', .1);
        end
    end
end

if length(iPatches)>1
    caxis([1 length(iPatches)])
end
colorbar;

view(30,40)
rotate3d;

set(gca,'Clipping','off');
end

