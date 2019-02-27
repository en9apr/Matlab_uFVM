function phiDiv = div(phi, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates the divergence of a vector/tensor field a
%--------------------------------------------------------------------------

% Check if vector
if size(phi(1, :, 1), 2)<3
    error('\n%s\n', 'The field values phi has to be of type vector/tensor');
end

numberOfElements = size(phi, 1);

% Determine if phi at an arbitrary element is a vector or a tensor
if size(phi(1, 1, :))<3
    % phi at each element is a vector
    phiDiv = zeros(numberOfElements, 1);
    
    % Calculate gradient of phi
    if nargin==1
        phiGrad = grad(phi);
    else
        phiGrad = varargin{1};
    end
    
    % Compute divergence of phi at each element
    for iElement=1:numberOfElements
        phiDiv(iElement) = phiGrad(iElement, 1, 1) + phiGrad(iElement, 2, 2) + phiGrad(iElement, 3, 3);
    end
else
    % phi at each element is a tensor
    phiDiv = zeros(numberOfElements, 3);
    
    % Calculate gradient of phi_x, phi_y and phi_z
    phiGrad1 = grad(phi(:, :, 1));
    phiGrad2 = grad(phi(:, :, 2));
    phiGrad3 = grad(phi(:, :, 3));
    for iElement=1:numberOfElements
        phiDiv(iElement, 1) =  phiGrad1(iElement, 1, 1) + phiGrad1(iElement, 1, 2) + phiGrad1(iElement, 1, 3);
        phiDiv(iElement, 2) =  phiGrad2(iElement, 2, 1) + phiGrad2(iElement, 2, 2) + phiGrad2(iElement, 2, 3);
        phiDiv(iElement, 3) =  phiGrad3(iElement, 3, 1) + phiGrad3(iElement, 3, 2) + phiGrad3(iElement, 3, 3);
    end
end


