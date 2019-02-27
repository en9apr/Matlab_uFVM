function cfdSetupMdotFields
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets up the Mdot field
%--------------------------------------------------------------------------
theMdotName = 'mdot_f';
cfdSetupMeshField(theMdotName,'Faces');