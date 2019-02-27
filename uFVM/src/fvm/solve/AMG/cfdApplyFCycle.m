function finalResidual = cfdApplyFCycle(gridLevel,smootherType,maxLevels,preSweep,postSweep,rrf)
%
%
%
while gridLevel<maxLevels
    cfdRestrict(gridLevel,smootherType,preSweep);
    gridLevel = gridLevel + 1;
end
%
% Smooth coarsest level
%
cfdSolveAlgebraicSystem(gridLevel,smootherType,postSweep);
%
%
cfdProlongate(gridLevel,smootherType,postSweep,rrf);
%
cfdRestrict(gridLevel-1,smootherType);
%
% Smooth coarsest level
%
cfdSolveAlgebraicSystem(gridLevel,smootherType,postSweep);
%
%
cfdProlongate(gridLevel,smootherType,postSweep,rrf);
%
cfdProlongate(gridLevel-1,smootherType,postSweep,rrf);
%
cfdRestrict(gridLevel-2,smootherType);
%
cfdRestrict(gridLevel-1,smootherType,preSweep);
%
% Smooth coarsest level
%
cfdSolveAlgebraicSystem(gridLevel,smootherType,postSweep);
%
%
cfdProlongate(gridLevel,smootherType,postSweep,rrf);
%
cfdProlongate(gridLevel-1,smootherType,postSweep,rrf);
%
cfdProlongate(gridLevel-2,smootherType,postSweep,rrf);
%
cfdRestrict(gridLevel-3,smootherType);
%
cfdRestrict(gridLevel-2,smootherType,preSweep);
%
cfdRestrict(gridLevel-1,smootherType,preSweep);
%
% Smooth coarsest level
%
cfdSolveAlgebraicSystem(gridLevel,smootherType,postSweep);
%
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