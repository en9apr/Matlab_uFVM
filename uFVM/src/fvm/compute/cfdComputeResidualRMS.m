function res = cfdComputeResidualRMS(residual)

res = 0;

numberOfElements = length(residual);
for iElement=1:numberOfElements
    localResidual = residual(iElement);
    res = res + localResidual^2;    
end
res = sqrt(res/numberOfElements);