function coefficients = cfdGetCoefficients(iLevel)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
if(nargin==0)
    iLevel = 1;
end   

global Domain;

coefficients =  Domain.coefficients{iLevel};


