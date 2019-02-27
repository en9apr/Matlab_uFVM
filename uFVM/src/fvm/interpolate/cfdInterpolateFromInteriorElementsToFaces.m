function phi_f = cfdInterpolateFromInteriorElementsToFaces(theInterpolationScheme,phi,mdotf)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
theSize = size(phi);
theNumberOfComponents = theSize(2);
if(theNumberOfComponents > 3)
    echo('**** Error *****');
    phi
    exit;
elseif(theNumberOfComponents==3)
    theType='Vector';
elseif(theNumberOfComponents==1)
    theType='Scalar';
end

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfFaces = theMesh.numberOfFaces;
numberOfBFaces = numberOfFaces-numberOfInteriorFaces;

interiorFaces = theMesh.faces(1:numberOfInteriorFaces);

iOwners = [interiorFaces.iOwner];
iNeighbours = [interiorFaces.iNeighbour];

if(strcmp(theInterpolationScheme,'Hyperbolic'))
    
    vol = [theMesh.elements.volume];
    for iComponent=1:theNumberOfComponents
        phi_f(:,iComponent) = (vol(iOwners)+vol(iNeighbours)).*phi(iOwners,iComponent).*phi(iNeighbours,iComponent)./(vol(iNeighbours).*phi(iOwners,iComponent)+vol(iOwners).*phi(iNeighbours,iComponent));
    end
elseif(strcmp(theInterpolationScheme,'Upwind'))
    
    vol = [theMesh.elements.volume];
    for iComponent=1:theNumberOfComponents
        phi_f(:,iComponent) = (vol(iOwners)+vol(iNeighbours)).*phi(iOwners,iComponent).*phi(iNeighbours,iComponent)./(vol(iNeighbours).*phi(iOwners,iComponent)+vol(iOwners).*phi(iNeighbours,iComponent));
    end
    
elseif(strcmp(theInterpolationScheme,'Average'))
    gf = [interiorFaces.gf]';
    for iComponent=1:theNumberOfComponents   
        phi_f(:,iComponent) = gf.*phi(iNeighbours,iComponent)+(1-gf).*phi(iOwners,iComponent);
    end

else 
    theInterpolationScheme
    exit;
end

%
% phi_f at boundary faces are simple set equal to the values of phi at boundary patch values
%
iBElement = 1;
for iBFace=numberOfInteriorFaces+1:numberOfFaces
    phi(numberOfElements+iBElement,:) = phi(theMesh.faces(iBFace).iOwner,:);
    iBElement = iBElement + 1;
end

for iComponent=1:theNumberOfComponents
    phi_f(numberOfInteriorFaces+1:numberOfFaces,iComponent) = phi(numberOfElements+1:numberOfElements+numberOfBFaces,iComponent);
end
