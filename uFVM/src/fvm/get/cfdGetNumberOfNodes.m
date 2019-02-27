function theNumberOfNodes = cfdGetNumberOfNodes;
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


mesh=cfdGetMesh;

theNumberOfNodes = length(mesh.nodes);