function iLevel = cfdAgglomerate(maxCoarseLevels)
%===================================================

%  written by the CFD Group @ AUB, Fall 2016
%===================================================

minNumberOfParents = 5;

iLevel = 1;
while iLevel<=maxCoarseLevels 
    iLevel = iLevel + 1;
    theNumberOfParents = cfdAgglomerateLevel(iLevel);
    assembleAgglomeratedLHS(iLevel);
    if theNumberOfParents<=minNumberOfParents
        break;
    end
end