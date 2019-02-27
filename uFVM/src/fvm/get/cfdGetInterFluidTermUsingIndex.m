function theInterFluidTerm = cfdGetInterFluidTermUsingIndex(iIFT)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global CFDEnv;


theInterFluidTerm=CFDEnv.interfluidTerms{iIFT};
