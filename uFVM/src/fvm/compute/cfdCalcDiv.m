function div = cfdCalcDiv(vec_f)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theMesh = cfdGetMesh;

div=zeros(1,numberOfElements);

for iFace=1:numberOfInteriorFaces
    Sf = theMesh.faces(iFace).Sf;
    iOwner = theMesh.faces(iFace).iOwner;
    iNeighbour = theMesh.faces(iFace).iNeighbour;

    div(iOwner) = div(iOwner) + vec_f(iFace)'*Sf;
    div(iNeighbour) = div(iNeighbour) - vec_f(iFace)'*Sf;
end


for iFace=numberOfInteriorFaces+1:numberOfFaces
    Sf = fvmFaces(iFace).Sf;
    iOwner = fvmFaces(iFace).iOwner;

    div(iOwner) = div(iOwner) + vec_f(iFace)'*Sf;
end

volumes = [theMesh.elements.volume];
div = div./volumes;
