function iComponent = cfdGetComponentFromName(theFieldName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


hasComponent = strfind(theFieldName,'(');
if(isempty(hasComponent))
    iComponent = 1;
else 
    iComponent = str2num(theFieldName(hasComponent+1));
end
