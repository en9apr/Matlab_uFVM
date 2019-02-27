function cfdPrintTermAttributes(theTerm)

if ~strcmp(theTerm.name, '0') && ~strcmp(theTerm.name, 'Source') && ~strcmp(theTerm.name, 'Pressure Gradient')
    fprintf('\n%s', ['Registering ', theTerm.name,' Term ...']);
    fprintf('\n%s', ['Expression: ', theTerm.expression]);
    fprintf('\n%s', ['Coefficient: ', theTerm.coefficientName]);
    fprintf('\n%s', ['Variable: ', theTerm.variableName]);
    fprintf('\n%s\n', ['Scheme: ', theTerm.scheme]);     
elseif strcmp(theTerm.name, 'Source')
    fprintf('\n%s', ['Registering ', theTerm.name,' Term ...']);
    fprintf('\n%s\n', ['Expression: ', theTerm.expression]);  
elseif strcmp(theTerm.name, 'Implicit Source')
    fprintf('\n%s', ['Registering ', theTerm.name,' Term ...']);
    fprintf('\n%s\n', ['Expression: ', theTerm.expression]);
elseif strcmp(theTerm.name, 'Pressure Gradient')
    fprintf('\n%s', ['Registering ', theTerm.name,' Term ...']);
    fprintf('\n%s\n', ['Expression: ', theTerm.expression]);    
end