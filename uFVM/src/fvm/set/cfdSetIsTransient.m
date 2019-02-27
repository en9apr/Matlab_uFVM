function cfdSetIsTransient(truefalse)

global Domain

if(nargin==0)
    truefalse = true;
end

Domain.isTransient = truefalse;

end