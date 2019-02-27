function cfdSetIsMoving (truefalse)

global Domain

if (nargin==0)
    truefalse= true;
end

Domain.isMoving= truefalse;
end