function grad_f = cfdInterpolateGradientsFromElementsToInteriorFaces(theInterpolationScheme,grad,phi,mdot_f)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

if nargin<3
    theNumberOfComponents = size(grad, 3);
else
    theNumberOfComponents = size(phi, 2);
end

theMesh = cfdGetMesh;

numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

iOwners = [theMesh.faces(1:numberOfInteriorFaces).iOwner]';
iNeighbours = [theMesh.faces(1:numberOfInteriorFaces).iNeighbour]';
gf = [theMesh.faces(1:numberOfInteriorFaces).gf]';

grad_f = zeros(numberOfInteriorFaces, 3, theNumberOfComponents);

if strcmp(theInterpolationScheme,'Average')
    for iComponent=1:theNumberOfComponents
        grad_f(:,1,iComponent) = (1-gf).*grad(iNeighbours,1,iComponent) + gf.*grad(iOwners,1,iComponent);
        grad_f(:,2,iComponent) = (1-gf).*grad(iNeighbours,2,iComponent) + gf.*grad(iOwners,2,iComponent);
        grad_f(:,3,iComponent) = (1-gf).*grad(iNeighbours,3,iComponent) + gf.*grad(iOwners,3,iComponent);
    end
    
elseif strcmp(theInterpolationScheme,'Upwind')
    pos = zeros(size(mdot_f));
    pos((mdot_f>0))=1;
    %
    grad_f(:,1) = pos.*grad(iNeighbours,1) + (1-pos).*grad(iOwners,1);
    grad_f(:,2) = pos.*grad(iNeighbours,2) + (1-pos).*grad(iOwners,2);
    grad_f(:,3) = pos.*grad(iNeighbours,3) + (1-pos).*grad(iOwners,3);
    
elseif strcmp(theInterpolationScheme,'Downwind')
    pos = zeros(size(mdot_f));
    pos((mdot_f>0))=1;
    %
    grad_f(:,1) = (1-pos).*grad(iNeighbours,1) + pos.*grad(iOwners,1);
    grad_f(:,2) = (1-pos).*grad(iNeighbours,2) + pos.*grad(iOwners,2);
    grad_f(:,3) = (1-pos).*grad(iNeighbours,3) + pos.*grad(iOwners,3);
    
elseif strcmp(theInterpolationScheme,'Average:Corrected') || strcmp(theInterpolationScheme,'Average:Rhie-Chow')
    for iComponent=1:theNumberOfComponents
        grad_f(:,1,iComponent) = (1-gf).*grad(iNeighbours,1,iComponent) + gf.*grad(iOwners,1,iComponent);
        grad_f(:,2,iComponent) = (1-gf).*grad(iNeighbours,2,iComponent) + gf.*grad(iOwners,2,iComponent);
        grad_f(:,3,iComponent) = (1-gf).*grad(iNeighbours,3,iComponent) + gf.*grad(iOwners,3,iComponent);
        
        d_CF = [theMesh.elements(iNeighbours).centroid]' - [theMesh.elements(iOwners).centroid]';
        dmag = cfdMagnitude(d_CF);
        
        e_CF(:,1) = d_CF(:,1)./dmag;
        e_CF(:,2) = d_CF(:,2)./dmag;
        e_CF(:,3) = d_CF(:,3)./dmag;
        
        local_grad_mag_f = (phi(iNeighbours,iComponent)-phi(iOwners,iComponent))./dmag;
        local_grad(:,1) = local_grad_mag_f.*e_CF(:,1);
        local_grad(:,2) = local_grad_mag_f.*e_CF(:,2);
        local_grad(:,3) = local_grad_mag_f.*e_CF(:,3);
        
        local_avg_grad_mag = dot(grad_f(:,:,iComponent)',e_CF')';
        local_avg_grad(:,1) = local_avg_grad_mag.*e_CF(:,1);
        local_avg_grad(:,2) = local_avg_grad_mag.*e_CF(:,2);
        local_avg_grad(:,3) = local_avg_grad_mag.*e_CF(:,3);
        
        grad_f(:,:,iComponent) = grad_f(:,:,iComponent) - local_avg_grad + local_grad;
    end
else
    theInterpolationScheme
    grad,
    phi,mdot
    exit;
end
