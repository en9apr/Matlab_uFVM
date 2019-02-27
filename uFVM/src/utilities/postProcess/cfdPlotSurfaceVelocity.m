function cfdPlotSurfaceVelocity(p,u1,u2)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%

figure;

axis off
axis equal

theMesh = cfdGetMesh;
theVelocityField = cfdGetMeshField('U');

colormap jet;

theNumberOfComponents = size(theVelocityField.phi,2);

for iComponent=1:theNumberOfComponents
    phi_n(:,iComponent) = cfdInterpolateFromElementsToNodes(theVelocityField.phi(:,iComponent));
end

theNumberOfNodes = theMesh.numberOfNodes;
theNumberOfFaces = theMesh.numberOfFaces;

for iNode=1:theNumberOfNodes
    v(iNode,:) = theMesh.nodes(iNode).centroid';
end

for iFace=1:theNumberOfFaces
    f{iFace,:} = theMesh.faces(iFace).iNodes;
end

f = cfdConvertCelltoMatrix(f);

plane = [p u1 u2];

cfdDrawMesh(v, f, 'FaceAlpha', 0, 'EdgeAlpha', .1);
cfdDrawPlane3d(plane, 'FaceAlpha', .5, 'EdgeAlpha', .5);

% compute the edge list
e = meshEdges(f);
edges = [v(e(:,1), :) v(e(:,2), :)];

% identify which edges cross the mesh
inds = isBelowPlane(v, plane);
edgeCrossInds = find(sum(inds(e), 2) == 1);

% compute one intersection point for each edge
intersectionPoints = intersectEdgePlane(edges(edgeCrossInds, :), plane);

for iIntersectionPt=1:size(intersectionPoints,1)
    iPoint1 = e(edgeCrossInds(iIntersectionPt),1);
    point1 = edges(edgeCrossInds(iIntersectionPt),1:3);
    
    iPoint2 = e(edgeCrossInds(iIntersectionPt),2);
    point2 = edges(edgeCrossInds(iIntersectionPt),4:6);
    
    intersectionPoint = intersectionPoints(iIntersectionPt,:);
    
    d1 = sqrt(sum((point1 - intersectionPoint).^2));
    d2 = sqrt(sum((intersectionPoint - point2).^2));
    
    gf = d2/(d1 + d2);
    
    for iComponent=1:theNumberOfComponents
        phi_p(iIntersectionPt,iComponent) = gf*phi_n(iPoint1,iComponent) + (1 - gf)*phi_n(iPoint2,iComponent);
    end
    
end

quiver3(intersectionPoints(:,1),intersectionPoints(:,2),intersectionPoints(:,3),phi_p(:,1),phi_p(:,2),phi_p(:,3),'b');

rotate3d;
view(30,40);

set(gca,'Clipping','off');