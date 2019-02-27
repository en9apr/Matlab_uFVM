function [initialResidual, finalResidual] = cfdSolveAlgebraicSystem(gridLevel,smootherType,numberOfIterations,rrf)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

if(nargin==3)
    rrf = 0.1;
end

theCoefficients = cfdGetCoefficients(gridLevel);
ac = theCoefficients.ac;
anb = theCoefficients.anb;
bc = theCoefficients.bc;
cconn = theCoefficients.cconn;
dphi = theCoefficients.dphi;
theNumberOfElements = theCoefficients.numberOfElements;
%
% Compute initial residual
%
residualsArray = cfdComputeResidualsArray(theCoefficients);
initialResidual = sum(abs(residualsArray))/theNumberOfElements;
finalResidual = initialResidual;
%
if strcmp(smootherType,'ILU')
    %
    % Factorize Ax=b (Apply incomplete upper lower decomposition)
    %
    [dc, rc] = cfdFactorizeILU(ac,anb,bc,cconn);
    %
    % Solve system for a number of times
    %
    iter = 1;
    while ((iter<=numberOfIterations) && (finalResidual>(rrf*initialResidual)))
        dphi = cfdSolveILU(ac,anb,bc,dc,rc,cconn,dphi);
        theCoefficients.dphi = dphi;
        
        residualsArray = cfdComputeResidualsArray(theCoefficients);
        finalResidual = sum(abs(residualsArray))/theNumberOfElements;
        %
        iter = iter + 1;
    end    
elseif strcmp(smootherType,'SOR')
    %
    % Compute initial residual
    %
    iter = 1;
    while ((iter<=numberOfIterations) && (finalResidual>(rrf*initialResidual)))
        dphi = cfdSolveSOR(ac,anb,bc,cconn,dphi);
        theCoefficients.dphi = dphi;
        %
        % Check if termination criterion satisfied
        %
        residualsArray = cfdComputeResidualsArray(theCoefficients);
        finalResidual = sum(abs(residualsArray))/theNumberOfElements;
        %
        iter = iter + 1;
     end
end

theCoefficients.dphi = dphi;
cfdSetCoefficients(theCoefficients,gridLevel);