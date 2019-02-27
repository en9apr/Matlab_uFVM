function phiDDt = DDt(phi, phi_old)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates the material derivative of a scalar/vector
%   field
%--------------------------------------------------------------------------

% Get time settings
isTransient = cfdIsTransient;
if isTransient
    deltaT = cfdGetDt;
    phiDdt = (phi - phi_old) / deltaT;
else
    phiDdt = zeros(size(phi));    
end

theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;

phiDDt = phiDdt + dot(U', grad(phi)')';


