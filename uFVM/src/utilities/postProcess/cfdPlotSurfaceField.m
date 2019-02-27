function cfdPlotSurfaceField(theFieldName,p,u1,u2,varargin)

theMesh = cfdGetMesh;

theFigure = figure;

axis off
axis equal

theFigure.Name = theFieldName;
colormap jet;

title(theFieldName);

if strcmp(theFieldName,'Ma')
    theMeshField = cfdGetMeshField('U');
    Vmag = cfdMagnitude(theMeshField.phi);
    
    theTemperatureField = cfdGetMeshField('T');
    T = theTemperatureField.phi;
    
    theGasConstantField = cfdGetMeshField('R');
    R = theGasConstantField.phi;
    
    theSpecificHeatRatioField = cfdGetMeshField('Cp');
    k = theSpecificHeatRatioField.phi;
    
    theMeshField.phi = Vmag./sqrt(k.*R.*T);
    theMeshField.type = 'Scalar';
else
    theMeshField = cfdGetMeshField(theFieldName);
end

theNumberOfComponents = size(theMeshField.phi,2);

for iComponent=1:theNumberOfComponents
    phi_n(:,iComponent) = cfdInterpolateFromElementsToNodes(theMeshField.phi(:,iComponent));
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

if nargin<5
    faceAlpha = 0.01;
else
    faceAlpha = varargin{2};
end
cfdDrawMesh(v, f, 'faceAlpha', faceAlpha, 'EdgeAlpha', .1);

% cfdDrawPlane3d(plane,'FaceAlpha', 0.1, 'EdgeAlpha', .3);

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

x = round(intersectionPoints(:,1), 9);
y = round(intersectionPoints(:,2), 9);
z = round(intersectionPoints(:,3), 9);

if theNumberOfComponents>1
    c = cfdMagnitude(phi_p);
else
    c = phi_p;
end

n = cross(u1, u2);

if isequal(n,[1 0 0]) || isequal(n,[-1 0 0])
    dt = delaunayTriangulation(y,z);
elseif isequal(n,[0 1 0]) || isequal(n,[0 -1 0])
    dt = delaunayTriangulation(x,z);
else
    dt = delaunayTriangulation(x,y);
end

h = patch('vertices',[x y z],'faces',dt.ConnectivityList,'FaceVertexCData',c,'FaceColor','interp','EdgeAlpha',0.1);

view(30,40);
rotate3d;

set(gca,'Clipping','off');

colorbar;


