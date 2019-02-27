function theConstantNames = cfdGetConstantNames;
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theConstantNames={};
theNumberOfConstants = length(Domain.constants);
for iConstant = 1:theNumberOfConstants
   theConstantNames{iConstant} = Domain.constants{iConstant}.name; 
end