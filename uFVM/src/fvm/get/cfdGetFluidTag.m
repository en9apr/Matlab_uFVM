function theFluidTag = cfdGetFluidTag(theFieldUserName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFieldName = cfdConvertName(theFieldUserName);

hasFluid = strfind(theFieldName,'_fluid');
if(isempty(hasFluid)) 
    theFluidTag = '';
else
    theFluidTag = theFieldName(hasFluid:hasFluid+7);
end
