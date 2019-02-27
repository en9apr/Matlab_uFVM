function cfdUpdateScale(theEquationUserName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function will compute the scale of an equation and store it in the
%   result in the equation field, for Vector equations one scale will be
%   computed for all components
%--------------------------------------------------------------------------

theField = cfdGetMeshField(theEquationUserName);

phi = theField.phi;

phiMax = max(max(phi));
phiMin = min(min(phi));
phiScale = max(max(phiMax - phiMin, abs(phiMax)));
if(phiScale < 1e-6)
    phiScale = 1;
end
    
theField.phiScale = phiScale;
theField.maximum = phiMax;
theField.minimum = phiMin;

cfdSetMeshField(theField);
