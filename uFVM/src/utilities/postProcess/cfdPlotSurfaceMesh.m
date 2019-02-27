function cfdPlotSurfaceMesh
figure;
%
%
axis off
axis equal
%
cameratoolbar(gcf,'Show'); % show the camera toolbar

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

h = patch('vertices', v, 'faces', f, 'FaceColor', 'red');

rotate3d;
view(30,40)
end