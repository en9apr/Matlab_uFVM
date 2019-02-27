function theSpecieName = cfdGetSpecieName(theName)


cPos = strfind(theName,':');
if(isempty(cPos))
    theSpecieName = '';
    
elseif(length(cPos)==1)
    
    theSpecieName = '';
    
elseif(length(cPos)==2)
    
    theSpecieName = theName(cPos(1)+1:cPos(2)-1);
    
else
    disp(['Error in converting name ->' theName]);
    stop;
end