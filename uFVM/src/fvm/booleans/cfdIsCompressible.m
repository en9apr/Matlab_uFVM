function status = cfdIsCompressible
%===================================================
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
% If at least one fluid is compressible then the simulation is assumed
% compressible as overall
%
status = false;
theFluidNames = cfdGetFluidNames;
for iFluid=1:length(theFluidNames)
    theFluidName = theFluidNames{iFluid};
    theFluid = cfdGetFluid(theFluidName);
    if strcmp(theFluid.mode, 'compressible')
        status = true;
    end
end