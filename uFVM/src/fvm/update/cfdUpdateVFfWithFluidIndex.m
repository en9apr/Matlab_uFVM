function cfdUpdateVFfWithFluidIndex(iFluid)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
% get mesh info
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements= theMesh.elements;
% ////////////////////////////////////////
% get needed fields
theFluidTag = cfdGetFluidTagUsingIndex(iFluid);
%
% get MdotField
%
theMdotField = cfdGetModel(['Mdot' theFluidTag]);
mdot_f=[];

if(isempty(theMdotField))
    theNumberOfFaces = cfdGetNumberOfAllFaces;
   mdot_f = ones(theNumberOfFaces,1);
else
  mdot_f = theMdotField.phi;  
end



% ////////////////////////////////////////

% interpolate fields to faces

% Get VF Field
%
theVFField = cfdGetModel(['VF' theFluidTag]);
VF=theVFField.phi;

theVFfField = cfdGetModel(['VF_f' theFluidTag]);
%
% Initialize VF_f at interior faces
%
interiorFacesIndices = find([fvmFaces(:).patchIndex]==0);
for iFace = interiorFacesIndices
    theFace = fvmFaces(iFace);
    iElement1 = theFace.element1;
    iElement2 = theFace.element2;
    theElement1 = fvmElements(iElement1);
    theElement2 = fvmElements(iElement2);
    if (mdot_f(iFace)>0)
        VF_f(iFace)=VF(iElement1);
    else
        VF_f(iFace)=VF(iElement2);
    end
end
theVFfField.phi = VF_f;
%
% Initialize VF_f_f at boundary faces
%
theNumberOfPatches = theMesh.numberOfPatches;
for iPatch=1:theNumberOfPatches
   theNumberOfPatchElements = theMesh.patchSizes(iPatch);
   % loop over patchFaces
    patchFacesIndices = find([fvmFaces(:).patchIndex] == iPatch);
    %
    VFb = theVFField.phiPatches{iPatch};
      
     %
    for k=1:theNumberOfPatchElements
      VF_b(k) = VFb(k);
%       ////////////////////////////////////////
      iBFace = patchFacesIndices(k);
      theVFfField.phi(iBFace) = VF_b(k);
%       /////////////////////////////////////////
    end
%     theVFfField.phiPatches{iPatch} = VF_b;
      
%
end

cfdSetModel(theVFfField);

end
