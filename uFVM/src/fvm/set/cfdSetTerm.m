function theEquation = cfdSetTerm(theEquationUserName,theTermName,varargin)
%
%
%
theNumberOfSubTerms = length(varargin);

theEquationName = cfdConvertName(theEquationUserName);
theEquation = cfdGetModel(theEquationName);
theTerm = cfdSetupTerm(theTermName);


theTerm.coefficient='';
theTerm.variable=theEquationUserName;
theTerm.scheme='DEFAULT';



for iSubTerm = 1:2:theNumberOfSubTerms
   theTerm =  setfield(theTerm,varargin{iSubTerm},varargin{iSubTerm+1});
end


if(strcmp(theTermName,'Diffusion'))
   theEquation.gamma = theTerm.coefficient; 
end
if(strcmp(theTermName,'Transient'))
   theEquation.rho = theTerm.coefficient; 
end
if(strcmp(theTermName,'Convection'))
   theEquation.rho = theTerm.coefficient; 
end

theNumberOfTerms = length(theEquation.terms);
theEquation.terms{theNumberOfTerms+1} =  theTerm ;

cfdSetModel(theEquation);