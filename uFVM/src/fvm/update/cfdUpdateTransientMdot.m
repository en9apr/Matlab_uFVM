function cfdUpdateTransientMdot(theMdotName)

theMdotMeshField = cfdGetMeshField(theMdotName,'Faces','Step0');
if ~isempty(theMdotMeshField)
    cfdSetMeshField(theMdotMeshField,'Step1')
end

