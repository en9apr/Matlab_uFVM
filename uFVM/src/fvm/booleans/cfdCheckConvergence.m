function isConverged = cfdCheckConvergence(convergenceCriterion)

isConverged = true;
theEquationNames = cfdGetEquationNames;
for iEquation=1:length(theEquationNames)
    tolerance = cfdGetResidualTolerance(theEquationNames{iEquation});
    if length(convergenceCriterion{iEquation})==1
        if convergenceCriterion{iEquation}>tolerance
            isConverged = false;
            return;
        end
    else
        for iComponentEquation=1:length(convergenceCriterion{iEquation})
            if convergenceCriterion{iEquation}(iComponentEquation)>tolerance
                isConverged = false;
                return;
            end
        end
    end
end