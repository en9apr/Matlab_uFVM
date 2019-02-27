function theEquationTerms = cfdReadEquationTerms(theFieldName, theEquation)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads equation terms and stores them
%--------------------------------------------------------------------------

global Domain;

% remove spaces. better manipulation
theEquation = strrep(theEquation, ' ', '');

% Add equality operator
if isempty(strfind(theEquation, '='))
    theEquation = strjoin({theEquation, '=0'});
end

% Add 0 term once required
if strcmp(theEquation(1), '+') || strcmp(theEquation(1), '-')
    theEquation = ['0', theEquation];
else
    theEquation = ['0+', theEquation];
end

if strcmp(theEquation(strfind(theEquation, '=')+1), '+') || strcmp(theEquation(strfind(theEquation, '=')+1), '-')
    str = strsplit(theEquation, '=');
    theEquation = [str{1}, '=0', str{2}];
else
    str = strsplit(theEquation, '=');
    theEquation = [str{1}, '=0+', str{2}];
end

% Split terms based on +, - and =
[C, matches_ic] = strsplit(theEquation, {'+', '-', '='});
iArithmeticOperator = 1;
for iM=1:length(matches_ic)
    if length(matches_ic{iM})>1
        for i=1:length(matches_ic{iM})
            matches{iArithmeticOperator} = matches_ic{iM}(i);
            iArithmeticOperator = iArithmeticOperator + 1;
        end
    else
        matches{iArithmeticOperator} = matches_ic{iM};
        iArithmeticOperator = iArithmeticOperator + 1;
    end
end

% Check if splitting didn't happen within a term
iTerm = 1;
for iC=1:length(C)
    if isempty(C{iC})
        continue;
    end
    i = 1;
    numberOfOP = length(strfind(C{iC}, '('));
    numberOfCP = length(strfind(C{iC}, ')'));
    while numberOfOP~=numberOfCP
        C{iC} = strjoin({C{iC}, C{iC+i}}, matches{iC+i-1});
        C{iC+i} = '';
        matches{iC+i-1} = '';
        numberOfOP = length(strfind(C{iC}, '('));
        numberOfCP = length(strfind(C{iC}, ')'));
        i = i + 1;
    end
    iTerm = iTerm + 1;
end

% Store in cell
iTerm = 1;
for iC=1:length(C)
    if isempty(C{iC}) || strcmp(C{iC}, '0')
        continue;
    end
    terms{iTerm} = C{iC};
    iTerm = iTerm + 1;
end

iArithmeticOperator = 1;
for iO=1:length(matches)
    if isempty(matches{iO}) || strcmp(matches{iO}, '=')
        continue;
    end
    arithmeticOperators{iArithmeticOperator} = matches{iO};
    iArithmeticOperator = iArithmeticOperator + 1;
end

% Get index of equality operator
equalOperatorIndex = strfind(theEquation, '=');

% Store in cell
iTerm_stored = 1;
for iTerm=1:length(terms)
    theEquationTerms{iTerm_stored}.name = '';
    theEquationTerms{iTerm_stored}.expression = terms{iTerm};
    if strcmp(arithmeticOperators{iTerm}, '+')
        if strfind(theEquation, terms{iTerm})<equalOperatorIndex
            theEquationTerms{iTerm_stored}.sign = 1;
        else
            theEquationTerms{iTerm_stored}.sign = -1;
        end
    elseif strcmp(arithmeticOperators{iTerm}, '-')
        if strfind(theEquation, terms{iTerm})<equalOperatorIndex
            theEquationTerms{iTerm_stored}.sign = -1;
        else
            theEquationTerms{iTerm_stored}.sign = 1;
        end
    end
    iTerm_stored = iTerm_stored + 1;
end


% Assign Term Types. Treated implicitly all the terms that include any of
% the operators ddt, div and laplacian applied in the following way
% assuming the equation is a general scalar equation phi:
% ddt(rho*phi) -> rho is an volume scalar field
% div(psi*phi) -> psi is a surface scalar field
% laplacian(gamma*phi) -> gamma is a surface scalar field
% All other terms are treated explicitly

for iTerm=1:length(theEquationTerms)
    % Check if implicit or Explicit
    expression = theEquationTerms{iTerm}.expression;
    attributes = getTermAttributes(expression, theFieldName);
        
    if strcmp(attributes.fvOptions, 'Implicit')        
        if ~isfield(attributes, 'operators')
            theEquationTerms{iTerm}.name = 'Implicit Source';
            continue;
        end
        if strcmp(attributes.operators{1}.name, 'ddt')
            theEquationTerms{iTerm}.name = 'Transient';
            theEquationTerms{iTerm}.operator = 'ddt';
            theEquationTerms{iTerm}.variableName = attributes.variable;
            theEquationTerms{iTerm}.coefficientName = attributes.coefficient;
            
            % Get the term scheme and store it
            ddtSchemes = Domain.foam.fvSchemes.ddtSchemes;
            theEquationTerms{iTerm}.scheme = '';
            for iDdtScheme=1:size(ddtSchemes, 1)
                if strcmp(strrep(theEquationTerms{iTerm}.expression, ' ', ''), strrep(ddtSchemes{iDdtScheme, 1}, ' ', ''))
                    theEquationTerms{iTerm}.scheme = ddtSchemes{iDdtScheme, 2};
                    break;
                end
            end
            if isempty(theEquationTerms{iTerm}.scheme)
                theEquationTerms{iTerm}.scheme = ddtSchemes{1, 2};
            end
        elseif strcmp(attributes.operators{1}.name, 'div')
            theEquationTerms{iTerm}.name = 'Convection';
            theEquationTerms{iTerm}.operator = 'div';
            theEquationTerms{iTerm}.variableName = attributes.variable;
            theEquationTerms{iTerm}.coefficientName = attributes.coefficient;
            
            % Get the term scheme and store it
            divSchemes = Domain.foam.fvSchemes.divSchemes;
            theEquationTerms{iTerm}.scheme = '';
            for iDivScheme=1:size(divSchemes, 1)
                if strcmp(strrep(theEquationTerms{iTerm}.expression, ' ', ''), strrep(divSchemes{iDivScheme, 1}, ' ', ''))
                    if strcmp(divSchemes{iDivScheme, 2}, 'Gauss upwind')
                        scheme = 'UPWIND';
                    elseif strcmp(divSchemes{iDivScheme, 2}, 'Gauss linear')
                        scheme = 'SOU';
                    else
                        scheme = 'SOU';
                    end
                    theEquationTerms{iTerm}.scheme = scheme;
                    break;
                end
            end
            if isempty(theEquationTerms{iTerm}.scheme)
                if strcmp(divSchemes{iDivScheme, 2}, 'Gauss upwind')
                    scheme = 'UPWIND';
                elseif strcmp(divSchemes{iDivScheme, 2}, 'Gauss linear')
                    scheme = 'SOU';
                else
                    scheme = 'SOU';
                end
                theEquationTerms{iTerm}.scheme = scheme;
            end
        elseif strcmp(attributes.operators{1}.name, 'laplacian')
            theEquationTerms{iTerm}.operator = 'laplacian';
            theEquationTerms{iTerm}.variableName = attributes.variable;
            theEquationTerms{iTerm}.coefficientName = attributes.coefficient;
            
%             if strcmp(theEquationTerms{iTerm}.variableName, 'U')
%                 theEquationTerms{iTerm}.name = 'Stress';
%             else
%                 theEquationTerms{iTerm}.name = 'Diffusion';
%             end
            theEquationTerms{iTerm}.name = 'Diffusion';
            
            % Get the term scheme and store it
            laplacianSchemes = Domain.foam.fvSchemes.laplacianSchemes;
            theEquationTerms{iTerm}.scheme = '';
            for iLaplacianScheme=1:size(laplacianSchemes, 1)
                if strcmp(strrep(theEquationTerms{iTerm}.expression, ' ', ''), strrep(laplacianSchemes{iLaplacianScheme, 1}, ' ', ''))
                    theEquationTerms{iTerm}.scheme = laplacianSchemes{iLaplacianScheme, 2};
                    break;
                end
            end
            if isempty(theEquationTerms{iTerm}.scheme)
                theEquationTerms{iTerm}.scheme = laplacianSchemes{1, 2};
            end
        end
    else
        theEquationTerms{iTerm}.name = 'Source';
        if ~isfield(attributes, 'operators')
            theEquationTerms{iTerm}.operators = '';
        else
            theEquationTerms{iTerm}.operators = attributes.operators;
        end
    end
end



