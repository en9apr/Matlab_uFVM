function isTransient = cfdIsTransient
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

if strcmp(Domain.time.type, 'Transient')
    isTransient = true;
else
    isTransient = false;
end