function theFluxes = cfdAssembleConvectionTermDCSTOIC(theEquationName,theTerm,theFluxes,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

if(nargin==3)
    iComponent= 1;
end

theMesh = cfdGetMesh;

numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
iFaces = 1:numberOfInteriorFaces;

theEquationMeshField = cfdGetMeshField(theEquationName);
phi = theEquationMeshField.phi;
phiGrad = theEquationMeshField.phiGradient;

theFluidTag = cfdGetFluidTag(theEquationName);
theMdotName = ['Mdot' theFluidTag];
theMdotField = cfdGetMeshField(theMdotName,'Faces');
mdot_f = theMdotField.phi(iFaces);


iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';
pos = zeros(size(mdot_f));
pos(mdot_f>0)=1;


% find indices of U and D cells
iUpwind   =  pos   .*iOwners + (1-pos).*iNeighbours;
iDownwind = (1-pos).*iOwners +  pos   .*iNeighbours;

% find phi_C, phi_D and calculate phi_U
phi_C = phi(iUpwind,iComponent);
phi_D = phi(iDownwind,iComponent);
rCD   = [theMesh.elements(iDownwind).centroid]'-[theMesh.elements(iUpwind).centroid]';

phi_U = phi_D - 2*dot(phiGrad(iUpwind,:,iComponent)',rCD')';

SMALL= 1e-6;
% calculate phi_tildaC
nominator = phi_C-phi_U;
denominator = phi_D-phi_U;
divideLoc = find(~((denominator<SMALL) & (denominator>-SMALL)));
phi_tildaC = ones(size(phi_C));
phi_tildaC(divideLoc) = nominator(divideLoc)./denominator(divideLoc);

% get phi_tildaf from STOIC function
phi_tildaf = zeros(size(phi_tildaC));
% lower UPWIND section
phi_tildaf = phi_tildaf + (phi_tildaC <= 0)                      .*(phi_tildaC);
% intermediate section
phi_tildaf = phi_tildaf + (phi_tildaC > 0)   .*(phi_tildaC < 0.2).*(3*phi_tildaC);
% CDS section
phi_tildaf = phi_tildaf + (phi_tildaC >= 0.2).*(phi_tildaC < 0.5).*(0.5*phi_tildaC + 0.5);
% SMART section
phi_tildaf = phi_tildaf + (phi_tildaC >= 0.5).*(phi_tildaC < 5/6).*(0.75*phi_tildaC + 3/8);
% DOWNWIND section
phi_tildaf = phi_tildaf + (phi_tildaC >= 5/6).*(phi_tildaC < 1)  .*(ones(size(phi_tildaC)));
% upper UPWIND section
phi_tildaf = phi_tildaf + (phi_tildaC >= 1)                      .*(phi_tildaC);

% calculate phi_f
phi_f = phi_tildaf.*(phi_D-phi_U) + phi_U;

% calculate correction
corr = mdot_f .* (phi_f - phi_C);

% apply deferred correction
theFluxes.FLUXTf(iFaces) = theFluxes.FLUXTf(iFaces) + corr;
