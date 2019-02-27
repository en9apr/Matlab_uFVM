function phiLaplacian = laplacian(phi, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates the laplacian of a scalar/vector field phi
%--------------------------------------------------------------------------

% Get mesh info
if nargin==1
    phiGrad = grad(phi);
else
    phiGrad = varargin{1};
end

phiLaplacian = div(phiGrad);

