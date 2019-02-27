function cfdPlotMesh(varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function plots the mesh
%--------------------------------------------------------------------------

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

if(nargin>0)
    phiNodes = varargin{1};      
    h = patch('vertices', v, 'faces', f, 'FaceVertexCData',phiNodes,'FaceColor','interp', 'EdgeAlpha', .1);
else    
    phiNodes = zeros(length(v),1);
    h = patch('vertices', v, 'faces', f, 'FaceVertexCData', phiNodes,'FaceAlpha', 0.1, 'EdgeAlpha', .4);
end

rotate3d;
set(gca,'Clipping','off');
view(30,40)
end

