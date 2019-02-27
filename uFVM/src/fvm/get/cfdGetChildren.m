function children = cfdGetChildren(iLevel)

theCoefficients = cfdGetCoefficients(iLevel-1);
parents = theCoefficients.parents;

theCoarseLevelCoefficients = cfdGetCoefficients(iLevel);
theNumberOfElements = theCoarseLevelCoefficients.numberOfElements;
children = cell(theNumberOfElements,1);

for iElement=1:length(parents)
    iParent = parents(iElement);    
    children{iParent} = [children{iParent} iElement];    
end