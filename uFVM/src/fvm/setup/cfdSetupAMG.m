function cfdSetupAMG(theEquationName,varargin)

theConvertedEquationName = cfdConvertName(theEquationName);
theEquation = cfdGetModel(theConvertedEquationName);

theMultigridSolver.isActive = true;
theMultigridSolver.cycleType = 'V-Cycle';
theMultigridSolver.maxCycles = '30';
theMultigridSolver.termination = '0.1';
theMultigridSolver.maxCoarseLevels = '10';
theMultigridSolver.preSweep = '1';
theMultigridSolver.postSweep = '3';
theMultigridSolver.nCycles = '30';

theNumberOfSubTerms = length(varargin);
for iSubTerm = 1:2:theNumberOfSubTerms
    theMultigridSolver = setfield(theMultigridSolver,varargin{iSubTerm},varargin{iSubTerm+1});
end

theEquation.multigrid = theMultigridSolver;

cfdSetModel(theEquation);