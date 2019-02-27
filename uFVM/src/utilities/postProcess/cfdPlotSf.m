function cfdPlotSf(iFigure)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%

figure(iFigure)
theMesh = cfdGetMesh;
%
%
cfdPlotMesh;
%

theNumberOfFaces = theMesh.numberOfFaces;

theFaceCentroids = [theMesh.faces(:).centroid]';

x = theFaceCentroids(:,1);
y = theFaceCentroids(:,2);
z = theFaceCentroids(:,3);    

iFaces = 1:theNumberOfFaces;

Sf = [theMesh.faces(iFaces).Sf]';
Sfx = Sf(:,1);
Sfy = Sf(:,2);
Sfz = Sf(:,3);

hold on;
quiver3(x,y,z,Sfx,Sfy,Sfz,2);
axis equal;
colorbar;


