function cfdInitializeEquation(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes the equation
%--------------------------------------------------------------------------

theField = cfdGetModel(theFieldName);
theType = theField.type;
%
%---------------------------------------------------
% Update interior elements
%---------------------------------------------------
%
% Evaluate the Field
%
phi = cfdComputeFormulaAtLocale(theField.ic,'Elements',theType);

theMeshField = cfdGetMeshField(theFieldName);

theMeshField.phi = phi;

cfdSetMeshField(theMeshField);
%
%---------------------------------------------------
% Update Boundary Conditions
%---------------------------------------------------
%
cfdInitializeBoundaryConditions(theFieldName);
%---------------------------------------------------
% Update Gradients 
%---------------------------------------------------
cfdUpdateGradient(theFieldName);

