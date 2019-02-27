function cfdCorrectDensityField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function corrects the density field
%--------------------------------------------------------------------------

% Get mesh info
thePressureCorrectionField = cfdGetMeshField('p');
p = thePressureCorrectionField.phi;

theDensityField = cfdGetMeshField('rho');

theDrhodpField = cfdGetMeshField('C_rho');
C_rho = theDrhodpField.phi;

rho = C_rho .* p;
theDensityField.phi = rho;
cfdSetMeshField(theDensityField);
