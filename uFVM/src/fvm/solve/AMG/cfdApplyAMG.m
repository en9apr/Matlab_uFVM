function [initialResidual, finalResidual] = cfdApplyAMG(theEquation)
%===================================================

%  written by the CFD Group @ AUB, Fall 2016
%===================================================
%
% Read multigrid settings
%
cycleType = theEquation.multigrid.cycleType;
maxCycles = eval(theEquation.multigrid.maxCycles);
rrf = eval(theEquation.multigrid.termination);
maxCoarseLevels = eval(theEquation.multigrid.maxCoarseLevels);
preSweep = eval(theEquation.multigrid.preSweep);
postSweep = eval(theEquation.multigrid.postSweep);
smootherType = theEquation.theSolutionMethods.smootherType;
%
% Build course grids
%
maxLevels = cfdAgglomerate(maxCoarseLevels);
%
% Calculate initial residual
%
theCoefficients = cfdGetCoefficients;
residualsArray = cfdComputeResidualsArray(theCoefficients);
initialResidual = sum(abs(residualsArray));
finalResidual = initialResidual;
%
gridLevel = 1;
nCycle = 1;
if(strcmp(cycleType,'V-Cycle'))
    while ((nCycle<=maxCycles)&&(finalResidual>rrf*initialResidual))
        %
        % Apply V-Cycle
        %
        finalResidual = cfdApplyVCycle(gridLevel,smootherType,maxLevels,preSweep,postSweep,rrf);
        %
        nCycle = nCycle + 1;
    end
elseif(strcmp(cycleType,'F-Cycle'))
    while ((nCycle<=maxCycles)&&(finalResidual>rrf*initialResidual))
        %
        % Apply F-Cycle
        %
        finalResidual = cfdApplyFCycle(gridLevel,smootherType,maxLevels,preSweep,postSweep,rrf);
        %
        nCycle = nCycle + 1;
    end
elseif(strcmp(cycleType,'W-Cycle'))
    while ((nCycle<=maxCycles)&&(finalResidual>rrf*initialResidual))
        %
        % Apply W-Cycle
        %
        finalResidual = cfdApplyWCycle(gridLevel,smootherType,maxLevels,preSweep,postSweep,rrf);
        %
        nCycle = nCycle + 1;
    end
end

theEquation.multigrid.nCycles = num2str(nCycle);
cfdSetModel(theEquation);