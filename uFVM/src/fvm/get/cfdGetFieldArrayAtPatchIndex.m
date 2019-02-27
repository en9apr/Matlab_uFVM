function phiPatch = cfdGetFieldArrayAtPatchIndex(theFieldName,thePatchIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theFieldBaseName = cfdGetBaseName(theFieldName);
theFluidName = cfdGetFluidName(theFieldName);
theFluid = cfdGetFluid(theFluidName);
theFluidTag = theFluid.tag;

theField=cfdGetModel(theFieldName);

if(isempty(theField))
    if(strcmp(theFieldBaseName,'Uplus'))
        phiPatch = cfdComputeUplus(thePatchIndex,theFluidTag);
    elseif(strcmp(theFieldBaseName,'Yplus'))
        phiPatch = cfdComputeYplus(thePatchIndex,theFluidTag);
    elseif(strcmp(theFieldBaseName,'Tauwall'))
        phiPatch = cfdComputeTauwall(thePatchIndex,theFluidTag);
    else
       disp (['****** ERROR ***** name not recognized: ' theFieldName]);
       stop;
    end
    
    
else
    
   phiPatch=theField.phiPatches{thePatchIndex};
    
end