function cfdCorrectFinerLevelSolution(gridLevel)

theCoefficients = cfdGetCoefficients(gridLevel);
DPHI = theCoefficients.dphi;

theFinerLevelCoefficients = cfdGetCoefficients(gridLevel-1);
dphi = theFinerLevelCoefficients.dphi;
theParents = theFinerLevelCoefficients.parents;
theNumberOfFineElements = theFinerLevelCoefficients.numberOfElements;

for iFineElement=1:theNumberOfFineElements
    iParent = theParents(iFineElement);
    dphi(iFineElement) = dphi(iFineElement) + DPHI(iParent);
end

theFinerLevelCoefficients.dphi = dphi;
cfdSetCoefficients(theFinerLevelCoefficients,gridLevel-1);