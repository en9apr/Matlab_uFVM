function cfdRunCase
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function runs the problem by iterating until convergence the
%   different equations previously defined by the user
%--------------------------------------------------------------------------

% Call time settings
time = cfdGetTime;
startTime = time.startTime;
endTime = time.endTime;
deltaT = time.deltaT;

% Run case
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);

cfdPrintHeader;

if cfdIsTransient
    timeIter = 1;
    cumulativeIter = 1;
    for t=startTime:deltaT:endTime
        % Time settings
        currentTime = t + deltaT;
        cfdSetCurrentTime(currentTime);
        
        % Update previous time step fields
        cfdTransientUpdate;
        
        cfdPrintCurrentTime(currentTime);
        
        % Loop until convergence for the current time step
        for iter=1:20
            cfdPrintIteration(cumulativeIter);
            cfdPrintResidualsHeader;
            %
            cfdUpdateFields;
            %
            for iEquation=1:theNumberOfEquations
                % Assemble the current equation and correct it
                [rmsResidual, maxResidual, lsResBefore, lsResAfter] = cfdAssembleAndCorrectEquation(theEquationNames{iEquation});
                
                % Print the equation residuals
                cfdPrintResiduals(cfdGetBaseName(theEquationNames{iEquation}),rmsResidual,maxResidual,lsResBefore,lsResAfter);
                
                % If multigrid solver is assigned, print the AMG solver
                % settings
                theEquation = cfdGetModel(theEquationNames{iEquation});
                isMultigrid = theEquation.multigrid.isActive;
                if isMultigrid
                    cfdPrintLinearSolver(theEquationNames{iEquation});
                    if iEquation<theNumberOfEquations
                        cfdPrintResidualsHeader;
                    end
                end
                
                % Store RMS residuals to check for convergence later on
                convergenceCriterion{iEquation}(1:length(maxResidual)) = rmsResidual;
                
                if any(isnan(rmsResidual))
                    msgbox('Divergence detected!');
                    return;
                end
            end
            fprintf('|==========================================================================|\n');
            cfdPrintCPUTime;
            
            cfdPlotRealTimeResiduals(cumulativeIter);
            
            % Check for convergence at each iteration at the current time
            % step
            isConverged = cfdCheckConvergence(convergenceCriterion);
            if isConverged
                fprintf('Solution is converged!\n')
                break;
            end
            cumulativeIter = cumulativeIter + 1;
        end
        cfdWriteResults(timeIter, currentTime);
        timeIter = timeIter + 1;
    end
    
else
    %
    % If steady state
    %
    currentIter = 1;
    for t=startTime:deltaT:endTime
        
        % Time settings
        currentTime = t + deltaT;
        cfdSetCurrentTime(currentTime);
        
        % Update transient for false transience update
        cfdTransientUpdate;
        
        cfdPrintIteration(currentIter);
        cfdPrintResidualsHeader;
        %
        cfdUpdateFields;
        %
        for iEquation=1:theNumberOfEquations
            % Assemble and Solve
            [rmsResidual, maxResidual, lsResBefore, lsResAfter] = cfdAssembleAndCorrectEquation(theEquationNames{iEquation});
            
            cfdPrintResiduals(cfdGetBaseName(theEquationNames{iEquation}),rmsResidual,maxResidual,lsResBefore,lsResAfter);
            
            theEquation = cfdGetModel(theEquationNames{iEquation});
            isMultigrid = theEquation.multigrid.isActive;
            if isMultigrid
                cfdPrintLinearSolver(theEquationNames{iEquation});
                if iEquation<theNumberOfEquations
                    cfdPrintResidualsHeader;
                end
            end
            
            convergenceCriterion{iEquation}(1:length(maxResidual)) = rmsResidual;
            
            if any(isnan(rmsResidual))
                msgbox('Divergence detected!');
                return;
            end
        end
        fprintf('|==========================================================================|\n');
        cfdPrintCPUTime;
        
        cfdWriteResults(currentIter, currentTime);
        
        cfdPlotRealTimeResiduals(currentIter);
        
        isConverged = cfdCheckConvergence(convergenceCriterion);
        if isConverged
            fprintf('Solution is converged!\n')
            break;
        end
        currentIter = currentIter + 1;
    end
end
