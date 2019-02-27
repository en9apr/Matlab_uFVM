function theFluidIndex = cfdGetFluidIndex(theFieldUserName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theFieldName = cfdConvertName(theFieldUserName);

hasFluid = strfind(theFieldName,'_fluid');
% if(isempty(hasFluid)) 
%     theFluidIndex = 0;
% else
    theFluidIndex = str2num(theFieldName(hasFluid+6:hasFluid+7));
% end
