function cfdUpdateGradientGauss0(theFieldName)
%===================================================
%
%  computes the gradient for a field at the centroids of 
%  the 'Elements'using a first order gauss interpolation
%  no correction for non-conjuntionality is applied
%  'TheFieldName' is the name of a field defined in the 
%  database
%
%  written by the CFD Group @ AUB, Fall 2013
%===================================================
%
theMeshField = cfdGetMeshField(theFieldName);
%
%
phiGrad = cfdComputeGradientGauss0(theMeshField.phi);


%-----------------------------------------------------
% Apply Under-relaxation to gradient 
%-----------------------------------------------------
%urfGrad = theField.urf*theField.urf;
%theMeshField.phiGradient = (1-urfGrad) *theMeshField.phiGradient + urfGrad*phiGrad;

theMeshField.phiGradient = phiGrad;

cfdSetMeshField(theMeshField);
