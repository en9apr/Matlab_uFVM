function theSpecieTag = cfdGetSpecieTag(theFieldUserName)

theFieldName = cfdConvertName(theFieldUserName);

hasSpecies = strfind(theFieldName,'_specie');
if(isempty(hasSpecies)) 
    theSpecieTag = '';
else
    theSpecieTag = theFieldName(hasSpecies:hasSpecies+7);
end
