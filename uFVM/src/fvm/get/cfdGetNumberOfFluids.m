function theNumberOfFluids = cfdGetNumberOfFluids
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

if(isempty(cfdGetFluids))
    theNumberOfFluids = 0;
else
    theNumberOfFluids = length(cfdGetFluids);
end