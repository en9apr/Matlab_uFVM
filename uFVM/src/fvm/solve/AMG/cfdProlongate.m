function cfdProlongate(gridLevel,smootherType,postSweep,rrf)
%
%
%
cfdCorrectFinerLevelSolution(gridLevel);
cfdSolveAlgebraicSystem(gridLevel-1,smootherType,postSweep,rrf);