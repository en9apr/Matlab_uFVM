function theFormula = modifyFormula(theFormula)

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
    phi_old = theMeshField_old.phi;
    
    argument_mod = [argument, ', phi_old'];
    theFormula = strrep(theFormula, argument, argument_mod);        
end
