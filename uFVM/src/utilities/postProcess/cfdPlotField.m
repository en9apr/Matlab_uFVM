function cfdPlotField(theFieldUserName, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function plots the contour of the field
%--------------------------------------------------------------------------

% Re-construct domain if not available
cfdReconstructDomain;

% Get the field name
theFieldName = cfdConvertFieldPhaseName(theFieldUserName);

if nargin>1
    theField = varargin{1};
else
    theField = cfdGetMeshField(theFieldName);
end

if strcmp(theField.type, 'Vector')
    phiNodes = cfdInterpolateFromElementsToNodes(sqrt(dot(theField.phi(:,:)',theField.phi(:,:)')));
    cfdPlotMesh(phiNodes);
    colorbar;
else
    phiNodes = cfdInterpolateFromElementsToNodes(theField.phi(:,1));
    cfdPlotMesh(phiNodes);
    colorbar;
end

title(theFieldUserName);
colormap jet