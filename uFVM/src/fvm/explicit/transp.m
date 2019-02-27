function transpA = transp(A)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function computes the transpose of a tensor field
%--------------------------------------------------------------------------

numberOfElements = size(A, 1);
transpA = zeros(size(A));
if size(A, 3)>1
    for iElement=1:numberOfElements
        for i=1:3
            for j=1:3
                transpA(iElement, i, j) = A(iElement, j, i);
            end
        end
    end
else
    error('\n%s\n', 'transp(A): A must be square matrix');
end