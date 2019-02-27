function S = cfdEvaluateNonstandardSourceTerm(theTerm, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function evaluates a nonstandard source term
%--------------------------------------------------------------------------

% Get mesh info
theMesh = cfdGetMesh;
iElements = 1:theMesh.numberOfElements;

% Loop over fields and evaluate them
theFormula = theTerm.formula;
theTermFields = identifyFields(theFormula);
for iTermField=1:length(theTermFields)
    theFieldName = theTermFields{iTermField};
    theMeshField = cfdGetMeshField(theFieldName);
    
    % If the required parameter is not a field, then it is a constant
    if isempty(theMeshField)
        theConstantValue = cfdGetConstant(theFieldName);
        tt = [theFieldName ' = theConstantValue;'];
        eval(tt);
    else
        % Execute the phi field
        tt = [theMeshField.name ' = theMeshField.phi;'];
        eval(tt);
    end
end

% Change arithmetic operators to be consistent for array operations
if isempty(strfind(theFormula, '.*')) && isempty(strfind(theFormula, './')) && isempty(strfind(theFormula, '.^'))
    theFormula = strrep(theFormula, '*', '.*');
    theFormula = strrep(theFormula, '/', './');
    theFormula = strrep(theFormula, '^', '.^');
end

% If ddt or DDt are included, then add the old values of the argument
if ~isempty(strfind(theFormula, 'DDt'))
    ddtOperatorIndex = strfind(theFormula, 'DDt');
    opIndecis = strfind(theFormula, '(');
    
    for iOpIndecis=1:length(opIndecis)
        if opIndecis(iOpIndecis)>ddtOperatorIndex
            opIndex = opIndecis(iOpIndecis);
            break;
        end
    end
    
    cpIndex = getClosingParanthesisIndex(theFormula, opIndex);
    argument = theFormula(opIndex+1:cpIndex-1);
    
    theMeshField_old = cfdGetMeshField(argument, 'Elements', 'Step1');
    if isempty(theMeshField_old)
        % If argument not available, then each of the fields of the
        % argument at older time step should be made available alone
        argument_old = argument;
        theTermFields = indentifyFields(argument);
        for iTermField=1:length(theTermFields)
            theFieldName = theTermFields{iTermField};
            theMeshField = cfdGetMeshField(theFieldName, 'Elements', 'Step1');
            tt = [theMeshField.name '_old = theMeshField.phi;'];
            eval(tt);
            argument_old = strrep(argument_old, theMeshField.name, [theMeshField.name '_old']);
        end
        argument_mod = [argument, ' ,', argument_old];
    else
        phi_old = theMeshField_old.phi;
        argument_mod = [argument, ', phi_old'];
    end
    theFormula = strrep(theFormula, argument, argument_mod);
end

S = eval(theFormula);

if nargin==2
    iComponent = varargin{1};
    S = S(iElements, iComponent);
else
    S = S(iElements, 1);
end