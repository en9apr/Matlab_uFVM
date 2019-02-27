function thePropertyNames = cfdGetPropertyNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

thePropertyNames = {};
theNumberOfFields = length(Domain.fields);
theNumberOfScalarFields = 0;
for iField=1:theNumberOfFields
    if(strcmp(Domain.fields{iField}.class,'Property'))
       theNumberOfScalarFields = theNumberOfScalarFields +1;
       thePropertyNames{theNumberOfScalarFields} = Domain.fields{iField}.name; 
    end
end