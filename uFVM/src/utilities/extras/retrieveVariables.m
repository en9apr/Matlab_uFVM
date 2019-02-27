function variables = retrieveVariables(expression)

C = textscan(expression, '%s', 'delimiter', {'+', '-', '*', '/', '(', ')'});

variables = {};
iVariable = 1;
for iC=1:length(C{1})
    if ~isempty(C{1}{iC}) && isnan(str2double(C{1}{iC}))        
        variables{iVariable} = strtrim(C{1}{iC});
        iVariable = iVariable + 1;
    end    
end