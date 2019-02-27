function theBaseFieldName = cfdRemoveComponentFromName(theFieldName)


hasComponent = strfind(theFieldName,'(');
if(isempty(hasComponent))
    theBaseFieldName = theFieldName;
else 
    theBaseFieldName = theFieldName(1:hasComponent-1);
end
