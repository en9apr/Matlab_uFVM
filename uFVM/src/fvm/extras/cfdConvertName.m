function theConvertedName = cfdConvertName(theName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

cPos = strfind(theName,':');
if(isempty(cPos))
    theConvertedName = theName;
    theFluidNames = cfdGetFluidNames;
    theNumberOfFluids = cfdGetNumberOfFluids;
    for iFluid=1:theNumberOfFluids
        theFluidName = theFluidNames{iFluid};
        theFluidTag = cfdGetFluidTagUsingIndex(iFluid);
        if(strcmp(theName,theFluidName))
            theConvertedName = theFluidTag;
        end
    end
    
elseif(length(cPos)==1)
    
    theBaseName = theName(1:cPos-1);
    theFluidName = theName(strfind(theName,':')+1:length(theName));
    theFluid = cfdGetFluid(theFluidName);
    theFluidTag = theFluid.tag;
    theConvertedName = [theBaseName theFluidTag];
    
elseif(length(cPos)==2)
    
    theBaseName = theName(1:cPos(1)-1);
    theSpecieName = theName(cPos(1)+1:cPos(2)-1);
    theSpecie = cfdGetSpecie(theSpecieName);
    theSpecieTag = theSpecie.tag;
    
    theFluidName = theName(cPos(2)+1:length(theName));
    theFluid = cfdGetFluid(theFluidName);
    theFluidTag = theFluid.tag;
    
    theConvertedName = [theBaseName theSpecieTag theFluidTag];
    
else
    disp(['Error in converting name ->' theName]);
    stop;
end