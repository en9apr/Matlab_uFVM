function attributes = getTermAttributes(term, equationName)


% Search for calculas operators
operators = {'ddt', 'div', 'laplacian', 'grad', 'transp', 'DDt'};

numberOfOperators = 1;
for iOperator=1:length(operators)
    operator = operators{iOperator};
    if strfind(term, operator)
        operatorIndecis = strfind(term, operator);
        for iOperatorIndecis=1:length(operatorIndecis)
            str = term;
            str(1:operatorIndecis(iOperatorIndecis)+length(operator)-1) = 'o';
            opIndecis = strfind(str, '(');
            opIndex = opIndecis(1);
            cpIndex = getClosingParanthesisIndex(term, opIndex);
            
            % Extract arguments
            argument = term(opIndex+1:cpIndex-1);
            
            attributes.operators{numberOfOperators}.name = operator;
            attributes.operators{numberOfOperators}.argument = argument;
            
            numberOfOperators = numberOfOperators + 1;
        end
    end
end

% If no operators found. It is usually regarded as a source term, however,
% only if the source term includes linear components such that Q = a*phi,
% then, if a is negative (given that it is at the RHS), then it is treated
% implicitly to enhance diagonal dominance, otherwise it is treated
% explicitly
if numberOfOperators==1
    variableIndexInTerm = strfind(term, equationName);
    if length(variableIndexInTerm)==1
        if variableIndexInTerm==1
            if strcmp(term(variableIndexInTerm+length(equationName)), '*') || strcmp(term(variableIndexInTerm+length(equationName)), '/')
                attributes.fvOptions = 'Implicit';
                return;
            else
                attributes.fvOptions = 'Explicit';
                return;
            end
        elseif variableIndexInTerm==length(term)-length(equationName)+1
            if strcmp(term(variableIndexInTerm-1), '*')
                attributes.fvOptions = 'Implicit';
                return;
            else
                attributes.fvOptions = 'Explicit';
                return;
            end
        else
            attributes.fvOptions = 'Explicit';
            return;
        end
    else
        attributes.fvOptions = 'Explicit';
        return;
    end
end

% Check if this term is to be treated implicitly or explicitly
% in the finite volume method. This is done by checking if any
% operator exist in the argument
if length(attributes.operators)>1
    attributes.fvOptions = 'Explicit';
else
    if strcmp(attributes.operators{1}.name, 'grad')
        attributes.fvOptions = 'Explicit';
    else
        attributes.fvOptions = 'Implicit';
    end
end

% If implicit fvOptions, extract the variable name and the coefficient
if strcmp(attributes.fvOptions, 'Implicit')
    attributes.variable = equationName;
    argument = strrep(attributes.operators{1}.argument, ' ', '');
    
    if strcmp(equationName, argument)
        attributes.coefficient = '1';
    else
        % remove flapping arithmetic operators
        variableIndex = strfind(argument, equationName);
        if length(variableIndex)>1
            variableIndex = variableIndex(end);
        end
        
        if variableIndex==1
            argument(variableIndex:variableIndex+length(equationName)) = [];
        elseif strcmp(argument(variableIndex-1), '(')
            if strcmp(argument(variableIndex+length(equationName)), ')')
                if variableIndex==2
                    argument(1:+length(equationName)+3) = [];
                else
                    argument(variableIndex-2:variableIndex+length(equationName)) = [];
                end
            else
                argument(variableIndex:variableIndex+length(equationName)) = [];
            end
        else
            argument(variableIndex-1:variableIndex+length(equationName)-1) = [];
        end
        
        coefficient = argument;
        attributes.coefficient = coefficient;
    end
end

