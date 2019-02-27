function cfdUpdateEquation(theFieldName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%---------------------------------------------------
% update Boundary Conditions
%---------------------------------------------------
%
cfdUpdateBoundaryConditions(theFieldName);
%---------------------------------------------------
% update Gradients 
%---------------------------------------------------
cfdUpdateGradient(theFieldName);

