function cfdRestrict(gridLevel,smootherType,preSweep)
%
% Restrict level to coarser level
%
if(nargin==3)
    cfdSolveAlgebraicSystem(gridLevel,smootherType,preSweep);
end
%
theCoefficients = cfdGetCoefficients(gridLevel);
residual = cfdComputeResidualsArray(theCoefficients);
%
cfdUpdateRHS(gridLevel+1,residual);