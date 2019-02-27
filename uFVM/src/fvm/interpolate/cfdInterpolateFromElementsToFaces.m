function phi_f = cfdInterpolateFromElementsToFaces(theInterpolationScheme,phi,mdot_f)
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
    theType = 'Vector';
elseif(theNumberOfComponents==1)
    theType = 'Scalar';
end;

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfFaces = theMesh.numberOfFaces;
numberOfBFaces = numberOfFaces-numberOfInteriorFaces;

interiorFaces = theMesh.faces(1:numberOfInteriorFaces);

iOwners = [interiorFaces.iOwner];
iNeighbours = [interiorFaces.iNeighbour];

if(strcmp(theInterpolationScheme, 'Hyperbolic'))    
    vol = [theMesh.elements.volume];
    for iComponent=1:theNumberOfComponents
        phi_f(:,iComponent) = (vol(iOwners)+vol(iNeighbours)).*phi(iOwners,iComponent).*phi(iNeighbours,iComponent)./(vol(iNeighbours).*phi(iOwners,iComponent)+vol(iOwners).*phi(iNeighbours,iComponent));
    end
elseif(strcmp(theInterpolationScheme, 'Upwind'))
    pos = zeros(size(mdot_f));
    pos(mdot_f>0) = 1;
    for iComponent=1:theNumberOfComponents
        phi_f(:,iComponent) = phi(iOwners,iComponent).*pos + phi(iNeighbours,iComponent).*(1 - pos);
    end
elseif(strcmp(theInterpolationScheme, 'Average'))
    gC = [interiorFaces.gf]';
    for iComponent=1:theNumberOfComponents        
        phi_f(:,iComponent) = gC.*phi(iNeighbours,iComponent) + (1 - gC).*phi(iOwners,iComponent);
    end
else 
    theInterpolationScheme
    exit;
end

%
% phi_f at boundary faces are simple set equal to the values of phi at boundary patch values
%
for iComponent=1:theNumberOfComponents
    phi_f(numberOfInteriorFaces+1:numberOfFaces,iComponent) = phi(numberOfElements+1:numberOfElements+numberOfBFaces,iComponent);
end
