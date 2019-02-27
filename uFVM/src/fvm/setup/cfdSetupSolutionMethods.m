function cfdSetupSolutionMethods(theEquationName,varargin)

theConvertedEquationName = cfdConvertName(theEquationName);
theEquation = cfdGetModel(theConvertedEquationName);

theSolutionMethods.smootherType = 'ILU';

theNumberOfSubTerms = length(varargin);
for iSubTerm = 1:2:theNumberOfSubTerms
    theSolutionMethods = setfield(theSolutionMethods,varargin{iSubTerm},varargin{iSubTerm+1});
end

theEquation.theSolutionMethods = theSolutionMethods;

cfdSetModel(theEquation);