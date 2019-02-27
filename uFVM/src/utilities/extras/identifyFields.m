function theTermFields = identifyFields(theFormula)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function identified the fields included in the provided formula
%--------------------------------------------------------------------------

% remove operators
operators = {'ddt', 'div', 'laplacian', 'grad', 'transp','DDt', 'log', 'sin', 'cos', 'tan', 'cot'};
for iOperator=1:length(operators)
   theFormula = strrep(theFormula, operators{iOperator}, '');
end

theFieldNames = cfdGetFieldNames;
theConstantNames = cfdGetConstantNames;
for iName=1:length(theFieldNames)
    theNames{iName} = theFieldNames{iName};
end
for iName=1:length(theConstantNames)
    theNames{length(theFieldNames)+iName} = theConstantNames{iName};
end

theTermFields = {};

% Sort according to field name lengths
for iFieldName=1:length(theFieldNames)
    theFieldName = theFieldNames{iFieldName};
    fieldNameLength = length(theFieldName);  
    lengths(iFieldName, 1) = fieldNameLength;
    lengths(iFieldName, 2) = iFieldName;
end

lengths = sortrows(lengths);
for i=1:size(lengths, 1)
    theNames{i} = theFieldNames{lengths(i, 2)};
end

theNames = fliplr(theNames);

numberOfTermFields = 1;
for iName=1:length(theNames)
    fieldIndex = strfind(theFormula, theNames{iName});
    if fieldIndex
        theTermFields{numberOfTermFields} = theNames{iName};
        theFormula = strrep(theFormula, theNames{iName}, '');
        numberOfTermFields = numberOfTermFields + 1; 
    end        
end





