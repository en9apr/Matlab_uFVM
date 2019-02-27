function cfdUpdateRHS(gridLevel,residual)

theCoefficients = cfdGetCoefficients(gridLevel-1);
theParents = theCoefficients.parents;
theNumberOfElements = theCoefficients.numberOfElements;

theCoarseLevelCoefficients = cfdGetCoefficients(gridLevel);
theNumberOfCoarseElements = theCoarseLevelCoefficients.numberOfElements;

BC = zeros(theNumberOfCoarseElements,1);

for iFineElement=1:theNumberOfElements    
    iParent = theParents(iFineElement);
    BC(iParent) = BC(iParent) + residual(iFineElement);
end

theCoarseLevelCoefficients.bc = BC;

cfdSetCoefficients(theCoarseLevelCoefficients,gridLevel);