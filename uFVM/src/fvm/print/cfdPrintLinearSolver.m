function cfdPrintLinearSolver(theEquationName)

theEquation = cfdGetModel(theEquationName);
if(~isempty(theEquation))
    isMultigrid = theEquation.multigrid.isActive;
    if isMultigrid
        cycleType = theEquation.multigrid.cycleType;
        nCycles = theEquation.multigrid.nCycles;
        smootherType = theEquation.theSolutionMethods.smootherType;
        theNumberOfAMGLevels = cfdGetNumberOfAMGLevels;
        fprintf('|--------------------------------------------------------------------------|\n');
        fprintf('|   Solver   |    Cycle    |    Levels   |     Cycles      |    Smoother   |\n');
        fprintf('|--------------------------------------------------------------------------|\n');
        nDigits1 = length(nCycles);
        spaces1 = '       ';
        if nDigits1>1
            spaces1 = spaces1(1:end-nDigits1+1);
        end
        
        nDigits2 = length(num2str(theNumberOfAMGLevels));
        spaces2 = '     ';
        if nDigits2>1
            spaces2 = spaces2(1:end-nDigits2+1);
        end
        
        fprintf(['|    AMG     |   ',cycleType,'   |      ',num2str(theNumberOfAMGLevels),spaces2,' |       ',nCycles,spaces1,'  |      ',smootherType,'      |\n']);
        fprintf('|- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |\n');
    end    
end
