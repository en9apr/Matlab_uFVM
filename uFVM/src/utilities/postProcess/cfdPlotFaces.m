function cfdPlotFaces(iFaces)

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
patch('vertices', v, 'faces', f, 'FaceVertexCData', phiNodes,'FaceAlpha', 0.05, 'EdgeAlpha', .4);

for iFace=iFaces
    theFace = theMesh.faces(iFace);
    theNodes = theMesh.nodes(theFace.iNodes);
    local_phi = iFace;
    XYZ = [theNodes.centroid]';
    patch(XYZ(:,1),XYZ(:,2),XYZ(:,3),local_phi,'FaceAlpha', 1, 'EdgeAlpha', .1);
end

if length(iFaces)>1
    caxis([1 length(iFaces)])
end
colorbar;

view(30,40)
rotate3d;

set(gca,'Clipping','off');
end