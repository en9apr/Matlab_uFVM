function theFluidName = cfdGetFluidName(theUserName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theFluid = cfdGetFluid(theUserName);
if (isempty(theFluid))
    cPos = strfind(theUserName,':');
    if(isempty(cPos))

        cfPos = strfind(theUserName,'_fluid');
        if(isempty(cfPos)) 
          theFluidName = '';
        else
            theFluidTag = theUserName(cfPos:cfPos+7);
            theFluid = cfdGetFluidUsingTag(theFluidTag);
            theFluidName = theFluid.userName;
        end
    elseif(length(cPos)==1)

        theFluidName = theUserName(strfind(theUserName,':')+1:length(theUserName));

    elseif(length(cPos)==2)

        theFluidName = theUserName(cPos(2)+1:length(theUserName));
    else
        disp(['Error in getFluidName name ->' theUserName]);
        stop;
    end
else
    theFluidName = theFluid.name;
end