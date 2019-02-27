function isConstant = cfdIsConstant(theFullConstantName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global CFDEnv;

isConstant = false;

theConstantName = cfdConvertFormula(theFullConstantName);

theNumberOfConstants = length(CFDEnv.constants);
for iConstant=1:theNumberOfConstants
    if(strcmp(theConstantName,CFDEnv.constants{iConstant}.name))
        isConstant = true;
        return
    end
end
