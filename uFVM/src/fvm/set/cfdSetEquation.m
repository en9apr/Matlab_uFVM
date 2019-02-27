function cfdSetEquation(theEquation)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theEquationIndex = 0;
theNumberOfEquations = length(Domain.equations);
for iEquation=1:theNumberOfEquations
    if(strcmp(theEquation.name,Domain.equations.name))
        Domain.equations{iEquation} = theEquation;
        theEquationIndex = iEquation;
    end
end
if(theEquationIndex==0)
    Domain.equations{theNumberOfEquations+1} = theEquation;
end


end