function theBaseName = cfdGetBaseName(theFieldUserName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFieldName = cfdConvertName(theFieldUserName);

hasSpecies = strfind(theFieldName,'_specie');
hasFluid = strfind(theFieldName,'_fluid');
if(isempty(hasSpecies))
    if(isempty(hasFluid)) 
    	theBaseName = theFieldName;
    else
        theBaseName = theFieldName(1:hasFluid-1);
    end

else 
    theBaseName = theFieldName(1:hasSpecies-1);
end

end
