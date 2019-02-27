function theScalarField = cfdUpdateGradientLS0(theScalarField)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
theMesh = cfdGetMesh;
%
%=====================================================
% INTERIOR Gradients 
%=====================================================
theNumberOfElements = theMesh.numberOfElements;

phi = theScalarField.phi;
%
urfGrad = theScalarField.urf*theScalarField.urf;

phiGradPrev= [];
if(isfield(theScalarField,'phiGradient'))
    phiGradPrev = theScalarField.phiGradient;
else
    phiGradPrev = zeros(theNumberOfElements,2);
    urfGrad = 1.0;
end
phiGrad = zeros(theNumberOfElements,2);
%-----------------------------------------------------
% INTERIOR FACES contribution to gradient 
%-----------------------------------------------------
theElements = theMesh.elements;
%
for iElement=1:theNumberOfElements
   %
   theElementC = theElements(iElement);
   theNeighbourIndices = theElementC.neighbours;
   
   phiC = phi(iElement);
   
   A11=0;
   A12=0;
   A21=0;
   A22=0;
   B11=0;
   B22=0;
   for iNeighbour = theNeighbourIndices
       theElementN = theElements(iNeighbour);
       
       phiN = phi(iNeighbour);
       CN = theElementN.centroid - theElementC.centroid;
       
       A11 = A11 + CN(1)*CN(1);
       A22 = A22 + CN(2)*CN(2);
       A12 = A12 + CN(1)*CN(2);
       A21 = A21 + CN(2)*CN(1);
       B11 = B11 + CN(1)*(phiN-phiC);
       B22 = B22 + CN(2)*(phiN-phiC);
       
   end
   %
   A=[A11 A12;A21 A22];
   B=[B11;B22];
   myGradient = A\B;
   
   phiGrad(iElement,:) = myGradient';
   
end

%-----------------------------------------------------
% Get Average Gradient by dividing with element volume 
%-----------------------------------------------------
fvmElements = theMesh.elements;

for iElement =1:theNumberOfElements
   theElement = fvmElements(iElement);
   % under relax the final gradient
   phiGrad(iElement,:) = urfGrad* phiGrad(iElement,:) + (1-urfGrad)*phiGradPrev(iElement,:);
end


theScalarField.phiGradient = phiGrad;
