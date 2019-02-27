function theFieldNames = cfdGetFieldNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theNumberOfFields = length(Domain.fields);

for iField=1:theNumberOfFields
   theFieldNames{iField} = Domain.fields{iField}.name; 
end