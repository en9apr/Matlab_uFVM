function cfdSetupTermCoefficientField(theEquationName, fvmOperator)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets up the term coefficient field
%--------------------------------------------------------------------------

if strcmp(fvmOperator, 'ddt')
    theScalarFieldName = ['rho_',theEquationName,'eq'];
    cfdSetupMeshField(theScalarFieldName);
elseif strcmp(fvmOperator, 'div')
    theScalarFieldName = ['psi_',theEquationName,'eq'];
    cfdSetupMeshField(theScalarFieldName, 'Faces');
elseif strcmp(fvmOperator, 'laplacian')
    theScalarFieldName = ['gamma_',theEquationName,'eq'];
    cfdSetupMeshField(theScalarFieldName, 'Faces');
end