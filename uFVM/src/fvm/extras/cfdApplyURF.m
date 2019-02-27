function cfdApplyURF(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function applies implicit under-relaxation factor
%--------------------------------------------------------------------------

theEquation = cfdGetModel(theEquationName);

urf = theEquation.urf;

theCoefficients = cfdGetCoefficients;
theCoefficients.ac = theCoefficients.ac/urf;

cfdSetCoefficients(theCoefficients);