function theConstant = cfdSetupConstant(theName,theValue)

%theConstantName = cfdConvertName(theName);
global Domain

if(isfield(Domain,'constants')==0)
   Domain.constants = {}; 
end
theConstant.name = theName;
theConstant.userName = theName;
theConstant.value = theValue;

ll=length(Domain.constants);
Domain.constants{ll+1}=theConstant;