function truefalse = cfdGetIsMoving
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


global Domain;

truefalse = Domain.time.isMoving;

end

