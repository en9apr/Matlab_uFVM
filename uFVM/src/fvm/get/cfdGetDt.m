function deltaT = cfdGetDt
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


global Domain;

deltaT = Domain.time.deltaT;

end