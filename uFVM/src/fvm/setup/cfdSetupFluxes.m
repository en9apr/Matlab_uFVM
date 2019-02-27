function theFluxes = cfdSetupFluxes
%
%
%
%
%

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;

%numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
%coefficients.Upper = zeros(1,numberOfInteriorFaces);
%coefficients.Lower = zeros(1,numberOfInteriorFaces);
%coefficients.Diag = zeros(1,numberOfElements);
%coefficients.u = [theMesh.faces(1:numberOfInteriorFaces).iOwner];
%coefficients.l = [theMesh.faces(1:numberOfInteriorFaces).iNeighbour];

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;

ac = zeros(numberOfElements,1);
ac_old = zeros(numberOfElements,1);
bc = zeros(numberOfElements,1);

dphi = zeros(numberOfElements,1);

cconn = cell(numberOfElements,1);
csize = zeros(numberOfElements,1);
anb = cell(numberOfElements,1);

for iElement=1:numberOfElements
   connectivity = theMesh.elements(iElement).iNeighbours;
   cconn{iElement}=connectivity;
   csize(iElement) = length(connectivity);
   anb{iElement} = zeros(1,csize(iElement));
end
coefficients.ac=ac;
coefficients.ac_old=ac_old;
coefficients.bc=bc;
coefficients.anb=anb;

coefficients.dphi=dphi;

coefficients.csize=csize;
coefficients.cconn = cconn;

cfdSetCoefficients(coefficients);