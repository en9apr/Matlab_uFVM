function cfdSetCoefficients(theCoefficients,iLevel)

global Domain;

if(nargin==1)
    iLevel = 1;
end

Domain.coefficients{iLevel} = theCoefficients;