function cfdGenerateCoarseGrids(theEquation)

maxCoarseLevels = eval(theEquation.multigrid.maxCoarseLevels);
coarseningBy = eval(theEquation.multigrid.coarseningBy);

theCoefficients = cfdGetCoefficients;
theMesh = cfdGetMesh;

AC{1} = theCoefficients.ac;
ANB{1} = theCoefficients.anb;
CCONN{1} = theCoefficients.cconn;

theCoefficients.multigrid.AC{1} = AC{1};
theCoefficients.multigrid.ANB{1} = ANB{1};
theMesh.multigrid.CCONN{1} = CCONN{1};

gridLevel = 2;
while gridLevel<=maxCoarseLevels+1
    [children{gridLevel,1}, parents{gridLevel,1}, CCONN{gridLevel,1}] = cfdAgglomeration(ANB{gridLevel-1,1},CCONN{gridLevel-1,1},coarseningBy);
    [AC{gridLevel,1}, ANB{gridLevel,1}] = cfdGetCoarseLevelCoefficients(children{gridLevel},parents{gridLevel},CCONN{gridLevel},AC{gridLevel-1},ANB{gridLevel-1},CCONN{gridLevel-1});    
    
    numberOfCoarseElements = length(AC{gridLevel});
    
    if numberOfCoarseElements<=3
        break;
    end
    
    theMesh.multigrid.numberOfCoarseElements(gridLevel) = numberOfCoarseElements;
    theMesh.multigrid.children{gridLevel} = children{gridLevel};
    theMesh.multigrid.CCONN{gridLevel} = CCONN{gridLevel};
    theMesh.multigrid.parents{gridLevel} = parents{gridLevel};
    theCoefficients.multigrid.AC{gridLevel} = AC{gridLevel};    
    theCoefficients.multigrid.ANB{gridLevel} = ANB{gridLevel};
    gridLevel = gridLevel + 1;
end

theMesh.multigrid.numberOfGridLevels = gridLevel - 1;

cfdSetMesh(theMesh);
cfdSetCoefficients(theCoefficients);



