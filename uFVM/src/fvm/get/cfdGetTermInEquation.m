function theTerm = cfdGetTermInEquation(theEquationName,theTermName,varargin)


theEquation = cfdGetModel(theEquationName);
theNumberOfTerms = length(theEquation.terms);
theTerm = '';

for iTerm=1:theNumberOfTerms    
    if strcmp(theTermName,theEquation.terms{iTerm}.name)
        if nargin==2
            theTerm = theEquation.terms{iTerm};
        else
            theTermExpression = varargin{1}; 
            if strcmp(theTermExpression, theEquation.terms{iTerm}.formula)
                theTerm = theEquation.terms{iTerm};
            end
        end
    end
end