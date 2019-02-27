function finalResidual = cfdApplyVCycle(gridLevel,smootherType,maxLevels,preSweep,postSweep,rrf)
%
% Resitriction
%
while gridLevel<maxLevels
    cfdRestrict(gridLevel,smootherType,preSweep);
    gridLevel = gridLevel + 1;
end
%
% Smoothening the coarsest level
%
cfdSolveAlgebraicSystem(gridLevel,smootherType,postSweep);
%
% Prolongation
%
while gridLevel>1
    cfdProlongate(gridLevel,smootherType,postSweep,rrf);
    gridLevel = gridLevel - 1;
end
%
% Calculate final residual
%
theCoefficients = cfdGetCoefficients;
residualsArray = cfdComputeResidualsArray(theCoefficients);
finalResidual = sum(abs(residualsArray));