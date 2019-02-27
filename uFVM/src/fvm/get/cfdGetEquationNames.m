function theEquationNames = cfdGetEquationNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theEquationNames={};
theNumberOfFields = length(Domain.fields);
theNumberOfScalarFields = 0;
for iField=1:theNumberOfFields
    if(strcmp(Domain.fields{iField}.class,'Equation'))
       theNumberOfScalarFields = theNumberOfScalarFields +1;
       theEquationNames{theNumberOfScalarFields} = Domain.fields{iField}.name; 
    end
end