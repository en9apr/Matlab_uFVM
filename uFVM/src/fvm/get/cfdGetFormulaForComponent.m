function theComponentFormula = cfdGetFormulaForComponent(thFormula,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


semiColumn = strfind(thFormula,';');
if(isempty(semiColumn)) 
    theComponentFormula = thFormula;
else
    if(iComponent==1)
        theComponentFormula = thFormula(2:semiColumn(1)-1);
    elseif(iComponent==2)
        theComponentFormula = thFormula(semiColumn(1)+1:semiColumn(2)-1);
    elseif(iComponent==3)
        theComponentFormula = thFormula(semiColumn(2)+1:length(thFormula)-1);
    end
end
