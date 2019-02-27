function theFluxes = cfdAssembleMixtureMdotTerm(theEquationName,theTerm)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%     This function assembles pressure correction equation.
%
%     1. Incompressible
%     U_f = v_f.S_f = [v]_f.S_f -  [DVOL]_f.(P_grad_f -[P_grad]_f).S_f + [DVOL]_f.([B]_f -[[B]]_f).S_f + [Dt]_f (U_Old_f -[v_old]_f.S_f) + (1-URF)(U_f -[v]_f.S_f)
%                                        ----------
%                     I                      II        III                       IV      V                     VI       VII                   VIII   IX
%     where only the underline part is linearized, the rest is treated explicitly
%
%     Fluxf = (rho_f*U_f)
%
%     2. Compressible
%
%     assemble the term velocity_f(old) rho_f(new) - velocity_f(old) rho_f(old)
%                                 (A)                              (B)
%     the rho is first written as a function of P
%     rho_f(new) = rho_f(old) + d(rho)/d(P)*[P(new) - P(old)]
%                 = d(rho)/d(P)*P(new) +  rho_f(old) - d(rho)/d(P)*P(old)
%
%     Fluxf = (U_f)*(1/RT_f)*P_f
%     this is simply equivalent to saying
%     assemble terms X
%      (mdot_f/rho_f)*(?rho/?p) P'
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
theNumberOfFaces = theMesh.numberOfFaces;
theNumberOfElements = theMesh.numberOfElements;

% Face fluxes
theFluxes.FLUXC1f(1:theNumberOfFaces, 1) = 0;
theFluxes.FLUXC2f(1:theNumberOfFaces, 1) = 0;
theFluxes.FLUXVf(1:theNumberOfFaces, 1)  = 0;
theFluxes.FLUXTf(1:theNumberOfFaces, 1)  = 0;

% Volume fluxes
theFluxes.FLUXCE(1:theNumberOfElements, 1) = 0;
theFluxes.FLUXCEOLD(1:theNumberOfElements, 1) = 0;
theFluxes.FLUXTE(1:theNumberOfElements, 1) = 0;

theMdotField = cfdGetMeshField('mdot_f', 'Faces');
theMdotField.phi(1:theNumberOfFaces, 1) = 0;

theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    % Get volume fractions
    theVFField = cfdGetMeshField(['alpha', num2str(iFluid)]);
    vf_f = cfdInterpolateFromElementsToFaces('Average', theVFField.phi);
    
    % Get density
    theRhoField = cfdGetMeshField(['rho', num2str(iFluid)]);
    rho_f = cfdInterpolateFromElementsToFaces('Average', theRhoField.phi);
    
    % Get mdot_f
    thePhasicMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)], 'Faces');
    
    % Assemble mdot_f at interior faces
    theCurrentFluxes = cfdAssembleMdotTermInterior(iFluid);
    
    % Assemble mdot_f at boundary faces
    theCurrentFluxes = cfdAssembleMdotTermBoundary(theCurrentFluxes,theEquationName,theTerm, iFluid);

    % Fix Pressure if Needed
    theCurrentFluxes = cfdAssembleMdotTermFixPressure(theCurrentFluxes);
    
    % Create phasic fluxes to be used in correcting the phasic mdot_f
    theFluxes.phasic_FLUXC1f(:, iFluid) = vf_f .* theCurrentFluxes.FLUXC1f;
    theFluxes.phasic_FLUXC2f(:, iFluid) = vf_f .* theCurrentFluxes.FLUXC2f;
    theFluxes.phasic_FLUXVf(:, iFluid) = vf_f .* theCurrentFluxes.FLUXVf;
    theFluxes.phasic_FLUXTf(:, iFluid) = vf_f .* theCurrentFluxes.FLUXTf;
    
    % Normalize by density and multiply by volume fraction
    theFluxes.FLUXC1f = theFluxes.FLUXC1f + vf_f .* theCurrentFluxes.FLUXC1f ./ rho_f;
    theFluxes.FLUXC2f = theFluxes.FLUXC2f + vf_f .* theCurrentFluxes.FLUXC2f ./ rho_f;
    theFluxes.FLUXVf = theFluxes.FLUXVf + vf_f .* theCurrentFluxes.FLUXVf ./ rho_f;
    theFluxes.FLUXTf = theFluxes.FLUXTf + vf_f .* theCurrentFluxes.FLUXTf ./ rho_f;
    
    % Calculate mixture mdot_f
    theMdotField.phi = theMdotField.phi + vf_f .* thePhasicMdotField.phi;
end
cfdSetMeshField(theMdotField);

% Store the fluxes in order to use in correcting mdot_f term later on
cfdSetFluxes(theFluxes);

% Assmeble into global Matrix
cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);
cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);

end


%===================================================
% INTERIOR Assembly
%===================================================
function theFluxes = cfdAssembleMdotTermInterior(iFluid)

% Get mesh info
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theMesh.numberOfBFaces;

iElements = 1:theMesh.numberOfElements;

iFaces = 1:numberOfInteriorFaces;
iBFaces = numberOfInteriorFaces+1:numberOfInteriorFaces+numberOfBFaces;

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

iBOwners = [theMesh.faces(iBFaces).iOwner]';

vol = [theMesh.elements.volume]';

% Initialize local face and element fluxes at interior faces and elements
% repectively
FLUXCf(iFaces, 1) = 0;
FLUXFf(iFaces, 1) = 0;
FLUXVf(iFaces, 1) = 0;
FLUXTf(iFaces, 1) = 0;

FLUXCE(iElements, 1) = 0;
FLUXCEOLD(iElements, 1) = 0;
FLUXTE(iElements, 1) = 0;

% Get DU fields and interpolate to faces
theD1Field = cfdGetMeshField('DU1');
theD2Field = cfdGetMeshField('DU2');
theD3Field = cfdGetMeshField('DU3');

d1_bar_f = cfdInterpolateFromElementsToFaces('Average',theD1Field.phi);
d1_bar_f = d1_bar_f(iFaces);

d2_bar_f = cfdInterpolateFromElementsToFaces('Average',theD2Field.phi);
d2_bar_f = d2_bar_f(iFaces);

d3_bar_f = cfdInterpolateFromElementsToFaces('Average',theD3Field.phi);
d3_bar_f = d3_bar_f(iFaces);

% Get first computed mdot_f
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)], 'Faces');
mdot_f_previous = theMdotField.phi(iFaces);

% Get velocity field and velocity model
theVelocityField = cfdGetMeshField('U');
U_bar_f = cfdInterpolateFromElementsToFaces('Average',theVelocityField.phi);
U_bar_f = U_bar_f(iFaces,:);

theVelocityEquation = cfdGetModel('U');
mdot_f_URF = theVelocityEquation.urf;

% Get rho field and assign density at faces as the convected one
theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho = theDensityField.phi;
rho_f = cfdInterpolateFromElementsToFaces('Upwind', rho, mdot_f_previous);
rho_f = rho_f(iFaces);

% Get pressure field and interpolate to faces
thePressureField = cfdGetMeshField('p');
p_C = thePressureField.phi(iOwners);
p_F = thePressureField.phi(iNeighbours);
p_grad_bar_f = cfdInterpolateGradientsFromElementsToInteriorFaces('Average',thePressureField.phiGradient,thePressureField.phi);
p_grad_f =  cfdInterpolateGradientsFromElementsToInteriorFaces('Average:Corrected',thePressureField.phiGradient,thePressureField.phi);

Sf = [theMesh.faces(iFaces).Sf]';
CF = [theMesh.faces(iFaces).CN]';
eCF = [theMesh.faces(iFaces).eCN]';
%
% Assemble Coefficients
%
DU1_f = d1_bar_f;
DU2_f = d2_bar_f;
DU3_f = d3_bar_f;
%
% assemble term I
%     rho_f [v]_f.Sf
%
U_bar_f = dot(U_bar_f(:,:)',Sf(:,:)')';
FLUXVf = FLUXVf + rho_f.*U_bar_f;
%
% Assemble term II and linearize it
%      - rho_f ([DPVOL]_f.P_grad_f).Sf
%
DUSf = [DU1_f.*Sf(:,1),DU2_f.*Sf(:,2),DU3_f.*Sf(:,3)]; % S'f
eDUSf = [DUSf(:,1)./cfdMagnitude(DUSf),DUSf(:,2)./cfdMagnitude(DUSf),DUSf(:,3)./cfdMagnitude(DUSf)];

DUEf = [cfdMagnitude(DUSf).*eCF(:,1)./dot(eCF(:,:)',eDUSf(:,:)')',cfdMagnitude(DUSf).*eCF(:,2)./dot(eCF(:,:)',eDUSf(:,:)')',cfdMagnitude(DUSf).*eCF(:,3)./dot(eCF(:,:)',eDUSf(:,:)')'];
geoDiff = cfdMagnitude(DUEf)./cfdMagnitude(CF);

DUTf = DUSf - DUEf;

FLUXCf = FLUXCf + rho_f.*geoDiff;
FLUXFf = FLUXFf - rho_f.*geoDiff;
FLUXVf  = FLUXVf  - rho_f.*dot(p_grad_f(iFaces,:)',DUTf(:,:)')';
%
%  assemble term III
%    rho_f ([P_grad]_f.([DPVOL]_f.Sf))
%
FLUXVf = FLUXVf + rho_f.*dot(p_grad_bar_f(iFaces,:)',DUSf(:,:)')';
%
% assemble terms IV and V
%    rho_f [DBVOL]_f.([B]_f -[[B]]_f).S_f
%
bodyForces = cfdGetBodyForces;
if ~isempty(bodyForces)
    for iBodyForce=1:length(bodyForces)
        bodyForceFormula = bodyForces{iBodyForce};
        theBodyForceTerm = cfdGetTermInEquation('U', 'Source', bodyForceFormula);
        B_C = cfdEvaluateNonstandardSourceTerm(theBodyForceTerm);
        B_bar_f = cfdInterpolateFromInteriorElementsToFaces('Average', B_C);
        
        rC = [theMesh.elements(iOwners).centroid theMesh.elements(iBOwners).centroid]';
        rF = [theMesh.elements(iNeighbours).centroid theMesh.faces(iBFaces).centroid]';
        rCF = rF - rC;
        
        S = [theMesh.faces.Sf]';
        
        gC = [theMesh.faces.gf]';
        
        B_bar_bar_C = zeros(theMesh.numberOfElements, 3);
        for iElement=iElements
            iNeighbourFaces = theMesh.elements(iElement).iFaces;
            for iNeighbourFace=iNeighbourFaces
                for iComponent=1:3
                    B_bar_bar_C(iElement, iComponent) = B_bar_bar_C(iElement, iComponent) + gC(iNeighbourFace) * dot(B_bar_f(iNeighbourFace, :)', rCF(iNeighbourFace, :)) * S(iNeighbourFace, iComponent);
                end
            end
            B_bar_bar_C(iElement, :) = B_bar_bar_C(iElement, :) ./ vol(iElement);
        end
        B_bar_bar_bar_f = cfdInterpolateFromInteriorElementsToFaces('Average', B_bar_bar_C);
        FLUXVf = FLUXVf + rho_f.*dot((B_bar_bar_bar_f(iFaces,:) - B_bar_f(iFaces,:))',DUSf(:,:)')';
    end
end
%
% assemble terms VI and VII
%    [Dt]_f (U_Old_f -[v_old]_f.S_f)
%
isTransient = cfdIsTransient;
if isTransient
    theDensityFieldOld = cfdGetMeshField(['rho', num2str(iFluid)], 'Elements', 'Step1');
    rho_old = theDensityFieldOld.phi;
    rho_old_f = cfdInterpolateFromElementsToFaces('Upwind', rho_old, mdot_f_previous);
    rho_old_f = rho_old_f(iFaces);
    
    theVelocityFieldOld = cfdGetMeshField('U', 'Elements', 'Step1');
    
    U_old_bar_f = cfdInterpolateFromElementsToFaces('Average', theVelocityFieldOld.phi);
    U_old_bar_f = U_old_bar_f(iFaces,:);
    
    U_old_bar_f = (dot(U_old_bar_f(:,:)',Sf(:,:)'))';
    
    theMdotFieldOld = cfdGetMeshField(['mdot_f', num2str(iFluid)], 'Faces', 'Step1');
    mdot_old_f = theMdotFieldOld.phi(iFaces);
    
    theDUT1Field = cfdGetMeshField('DUT1');
    theDUT2Field = cfdGetMeshField('DUT2');
    theDUT3Field = cfdGetMeshField('DUT3');
    
    DUT = mean([theDUT1Field.phi,theDUT2Field.phi,theDUT3Field.phi],2);
    DUT_f = cfdInterpolateFromElementsToFaces('Average', DUT);
    DUT_f = DUT_f(iFaces);
    
    FLUXVf = FLUXVf + DUT_f.*(mdot_old_f - rho_old_f.*U_old_bar_f);      
end
%
% assemble terms VIII and IX
%     (1-URF)(U_f -[v]_f.S_f)
%
FLUXVf = FLUXVf + (1 - mdot_f_URF)*(mdot_f_previous - rho_f.*U_bar_f);
%
% compute Rhie-Chow interpolation of mdot_f and updated it in the data base
%
mdot_f = FLUXCf .* p_C + FLUXFf .* p_F + FLUXVf;
theMdotField.phi(iFaces) = mdot_f;
cfdSetMeshField(theMdotField);

FLUXTf = mdot_f;
%
% assemble terms X (for compressible flow)
%
applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'compressible')
    theDrhodpField = cfdGetMeshField('C_rho');
    C_rho = theDrhodpField.phi;
    C_rho_f = cfdInterpolateFromElementsToFaces('Upwind', C_rho, mdot_f);
    C_rho_f = C_rho_f(iFaces);
    
    FLUXCf = FLUXCf + (C_rho_f ./ rho_f) .* max(mdot_f, 0);
    FLUXFf = FLUXFf - (C_rho_f ./ rho_f) .* max(-mdot_f, 0);
    
    % Add transient contribution
    if isTransient
        deltaT = cfdGetDt;        
        
        FLUXCE = vol(iElements) .* C_rho(iElements) / deltaT;
        FLUXTE = (rho(iElements) - rho_old(iElements)) .* vol(iElements) / deltaT;
    end
end

%
% Assemble in Global Fluxes
%
theFluxes.FLUXC1f(iFaces,1) = FLUXCf;
theFluxes.FLUXC2f(iFaces,1) = FLUXFf;
theFluxes.FLUXVf(iFaces,1)  = FLUXVf;
theFluxes.FLUXTf(iFaces,1)  = FLUXTf;

theFluxes.FLUXCE(iElements,1) = FLUXCE;
theFluxes.FLUXCEOLD(iElements, 1) = FLUXCEOLD;
theFluxes.FLUXTE(iElements,1) = FLUXTE;

end



%===================================================
% FIX PRESSURE
%===================================================
function theFluxes = cfdAssembleMdotTermFixPressure(theFluxes)

theMesh = cfdGetMesh;
theCoefficients = cfdGetCoefficients;

iFixedElement = cfdGetFixedElement;
if(iFixedElement>0) && (iFixedElement<theMesh.numberOfElements)
    theMesh = cfdGetMesh;
    theElement = theMesh.elements(iFixedElement);
    iNeighbours = theElement.iNeighbours;
    
    % Set Neighbouring Coefficients to 0
    for iNBElement = iNeighbours
        theCoefficients.anb{iFixedElement,iNBElement} = 0;
    end
    theCoefficients.bc(iFixedElement) = 0;
end

% Update Cofficients
cfdSetCoefficients(theCoefficients);

end


function theFluxes = cfdAssembleMdotTermBoundary(theFluxes,theEquationName,theTerm, iFluid)

theEquation = cfdGetModel(theEquationName);
theMesh = cfdGetMesh;
iBFaces = theMesh.numberOfInteriorFaces+1:theMesh.numberOfFaces;

% Initialize all boundary fluxes
theFluxes.FLUXC1f(iBFaces,1) =  0;
theFluxes.FLUXC2f(iBFaces,1) =  0;
theFluxes.FLUXVf(iBFaces,1)  =  0;
theFluxes.FLUXTf(iBFaces,1)  =  0;

theNumberOfPatches = theMesh.numberOfPatches;
for iPatch=1:theNumberOfPatches
    % find the Physical Type
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    % find type of BC
    theBCType = theEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'noSlip')
            theFluxes =  cfdAssembleMdotTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        elseif strcmp(theBCType,'slip')
            theFluxes =  cfdAssembleMdotTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        elseif strcmp(theBCType,'zeroGradient')
            theFluxes =  cfdAssembleMdotTermWallZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        else
            error([theBCType '<<<< not implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            theFluxes =  cfdAssembleMdotTermInletInletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        elseif(strcmp(theBCType,'Supersonic Inlet'))
            theFluxes =  cfdAssembleMdotTermInletSupersonicInletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        elseif strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleMdotTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        else
            error([theBCType '<<<< not implemented']);
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'outlet')
            theFluxes =  cfdAssembleMdotTermOutletOutletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        elseif strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleMdotTermOutletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iFluid);
        else
            error([theBCType '<<<< not implemented']);
        end
        %
        % ERROR
        %
    elseif(strcmp(thePhysicalType,'empty'))
        
    elseif(strcmp(thePhysicalType,'symmetry'))
        
    else
        error([thePhysicalType '<<<< not implemented']);
    end
    %
end

end

%===================================================
% WALL- slip Conditions
%===================================================
function theFluxes =  cfdAssembleMdotTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0.0;
theFluxes.FLUXC2f(iBFaces) = 0.0;
theFluxes.FLUXVf(iBFaces)  = 0.0;
theFluxes.FLUXTf(iBFaces)  = 0.0;

end

%===================================================
% WALL-noslip Condition
%===================================================
function theFluxes =  cfdAssembleMdotTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0.0;
theFluxes.FLUXC2f(iBFaces) = 0.0;
theFluxes.FLUXVf(iBFaces)  = 0.0;
theFluxes.FLUXTf(iBFaces)  = 0.0;

end

%===================================================
% WALL-zeroGradient Condition
%===================================================
function theFluxes =  cfdAssembleMdotTermWallZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0.0;
theFluxes.FLUXC2f(iBFaces) = 0.0;
theFluxes.FLUXVf(iBFaces)  = 0.0;
theFluxes.FLUXTf(iBFaces)  = 0.0;

end

%===================================================
% INLET-Inlet Boundary Conditions
%===================================================
function theFluxes =  cfdAssembleMdotTermInletInletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

% Get mesh info
theMesh = cfdGetMesh;
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
%
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
%
%
Sf = [theMesh.faces(iBFaces).Sf]';
%
% Get Arrays
%
theVelocityField = cfdGetMeshField('U');
U_b = theVelocityField.phi(iBElements,:);

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi(iBElements);

thePressureField = cfdGetMeshField('p');
p_C = thePressureField.phi(iOwners);
p_b = thePressureField.phi(iBElements);

% ---------
%  STEP 1
% ---------
%  Assemble FLUXCb_1, FLUXCb_2 and FLUXVb coefficients
%
%---------------------------------------------------
% assemble RHIE-CHOW Interpolation Term I
%---------------------------------------------------
%
FLUXCb1(1:numberOfBFaces,1) =  0;
FLUXCb2(1:numberOfBFaces,1) =  0;
FLUXVb(1:numberOfBFaces,1)  =  0;
FLUXTb(1:numberOfBFaces,1)  =  0;
%
% ONLY assemble term I
%
U_b = dot(Sf(:,:)',U_b(:,:)')';
FLUXVb = FLUXVb + rho_b.*U_b;  % term I
%
% Store mdot_b
%
mdot_b = FLUXCb1.*p_C + FLUXCb2.*p_b + FLUXVb;
theMdotField.phi(iBFaces) = mdot_b;
cfdSetMeshField(theMdotField);

FLUXTb = mdot_b;
%
% assemble terms X (for compressible)
%
applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'compressible')
    theDrhodpField = cfdGetMeshField('C_rho');
    C_rho_b = theDrhodpField.phi(iBElements);
    
    FLUXCb1 = FLUXCb1 + (C_rho_b ./ rho_b) .* mdot_b;
end

theFluxes.FLUXC1f(iBFaces) =  FLUXCb1;
theFluxes.FLUXC2f(iBFaces) =  FLUXCb2;
theFluxes.FLUXVf(iBFaces)  =  FLUXVb;
theFluxes.FLUXTf(iBFaces)  =  FLUXTb;

end


%===================================================
% SUPERSONIC INLET-Inlet Boundary Conditions
%===================================================
function theFluxes =  cfdAssembleMdotTermInletSupersonicInletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)
%
theMesh = cfdGetMesh;
%
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
%
Sf = [theMesh.faces(iBFaces).Sf]';
%
% Get Arrays
%
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_previous_b = theMdotField.phi(iBFaces);
mdot_b = zeros(size(mdot_previous_b));

theVelocityField = cfdGetMeshField('U');
U_b = theVelocityField.phi(iBElements,:);

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi(iBElements);

thePressureField = cfdGetMeshField('p');
p = thePressureField.phi(iOwners);
p_b = thePressureField.phi(iBElements);

% ---------
%  STEP 1
% ---------
%  Assemble FLUXCb_1, FLUXCb_2 and FLUXVb coefficients
%
%---------------------------------------------------
% assemble RHIE-CHOW Interpolation Term I
%---------------------------------------------------												%
% loop over patchFaces
%
lBFaces = 1:theNumberOfPatchFaces;
%for iBFace=1:theNumberOfPatchFaces
%
local_FLUXCb1(lBFaces,1) =  0.0;
local_FLUXCb2(lBFaces,1) =  0.0;
local_FLUXVb(lBFaces,1)  =  0.0;
%
% ONLY assemble term I
%
U_b = dot(Sf(lBFaces)',U_b(lBFaces)')';
local_FLUXVb = local_FLUXVb + rho_b(lBFaces).*U_b;  % term I
%
% assemble terms X (For Compressible Flows)
%     (mdot_f/rho_f)*(?rho/?p) P'
% where mdot_f is the newly computed mdot_f
%
local_mdot_b = local_FLUXCb1.*p(lBFace) + local_FLUXCb2.*p_b(lBFace) + local_FLUXVb;
%
%local_FLUXCb1 = local_FLUXCb1 + 0.0;
%local_FLUXCb2 = local_FLUXCb2 + (local_mdot_b/rho_b[iBFace])*drhodp_b[iBFace];
%local_FLUXVb   = local_FLUXVb -(0.0*(pressure[iElement1]+ Pref) + (local_mdot_b/rho_b[iBFace])*drhodp_b[iBFace]*(pressure_b[iBFace]+ Pref));
% local_FLUXVb = local_FLUXVb -local_mdot_b;
%
theFluxes.FLUXC1f(iBFaces) = theFluxes.FLUXC1f(iBFaces) + local_FLUXCb1;
theFluxes.FLUXC2f(iBFaces) = theFluxes.FLUXC2f(iBFaces) + local_FLUXCb2;
theFluxes.FLUXVf(iBFaces)  = theFluxes.FLUXVf(iBFaces) + local_FLUXVb;
%
local_FLUXTb = local_mdot_b;
theFluxes.FLUXTf(iBFaces)  = local_FLUXTb;
mdot_b(lBFaces) = mdot_b(lBFaces) + local_FLUXTb;
theMdotField.phi(iBFaces) = mdot_b;
cfdSetMeshField(theMdotField);
end


%===================================================
% Inlet-fixedValue Condition
%===================================================
function theFluxes =  cfdAssembleMdotTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

%
theMesh = cfdGetMesh;
%
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
%
Sb = [theMesh.faces(iBFaces).Sf]';
CN = [theMesh.faces(iBFaces).CN]';
eCN = [theMesh.faces(iBFaces).eCN]';
%
% Get Arrays
%
theMdotName = ['mdot_f', num2str(iFluid)];
theMdotField = cfdGetMeshField(theMdotName,'Faces');
mdot_previous_b = theMdotField.phi(iBFaces);
mdot_b = zeros(size(mdot_previous_b));

theVelocityField = cfdGetMeshField('U');
U_b = theVelocityField.phi(iBElements,:);

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi(iBElements);

thePressureField = cfdGetMeshField('p');
pressure = thePressureField.phi(iOwners);
pressure_b = thePressureField.phi(iBElements);
pressure_grad = thePressureField.phiGradient(iOwners);
pressure_grad_b=thePressureField.phiGradient(iBElements);
%
theD1Field = cfdGetMeshField('DU1');
theD2Field = cfdGetMeshField('DU2');
theD3Field = cfdGetMeshField('DU3');
%
d1_b = theD1Field.phiPatches{iPatch};
d2_b = theD2Field.phiPatches{iPatch};
d3_b = theD3Field.phiPatches{iPatch};
%
% ---------
%  STEP 1
% ---------
%  Assemble the Coefficients FLUXCb_1, FLUXCb_2 and FLUXVb
%
%-------------------------------------------------------------
% assemble RHIE-CHOW Interpolation Terms I,II,III,VIII and IX
%-------------------------------------------------------------
%
% loop over patchFaces
%
for k=1:numberOfBFaces
    %
    %
    local_FLUXCb1 =  0.0;
    local_FLUXCb2 =  0.0;
    local_FLUXVb =  0.0;
    local_FLUXTb =  0.0;
    %
    % Assemble Coefficients
    %
    DP_b = vf_b(iBFace)*[d1_b(k),0.0,0.0;
        0.0,d2_b(k),0.0;
        0.0,0.0,d3_b(k)];
    DUE = [eCN(k,1).*cfdMagnitude(DP_b),eCN(k,2).*cfdMagnitude(DP_b),eCN(k,3).*cfdMagnitude(DP_b)];
    %
    % Assemble term I
    %
    U_bar_b = U_b(k)*Sb(k);
    local_FLUXVb =  local_FLUXVb + rho(k)*U_bar_b;
    %
    % Assemble term II and linearize it
    %
    geoDiff = (Sb(k,:)*(DP_b*Sb(k,:)'))/(CN(k,:)*Sb(k,:)');
    local_FLUXCb1 =  local_FLUXCb1 + rho_b(k)*geoDiff;
    local_FLUXCb2 = local_FLUXCb2 - rho_b(k)*geoDiff;
    local_FLUXVb  = local_FLUXVb - rho_b(k)*((pressure_grad_b(k,:)')'*DUE(k)');
    %
    %  assemble term III
    %
    local_FLUXVb = local_FLUXVb + rho_b(iBFace)*((DP_b*pressure_grad(k,:)')'*Sb(k)');
    %
    
    theFluxes.FLUXC1f(iBFaces(iBFace)) = theFluxes.FLUXC1f(iBFaces(iBFace))  + local_FLUXCb1;
    theFluxes.FLUXC2f(iBFaces(iBFace)) = theFluxes.FLUXC2f(iBFaces(iBFace))  + local_FLUXCb2;
    theFluxes.FLUXVf(iBFaces(iBFace)) = theFluxes.FLUXVf(iBFaces(iBFace))  + local_FLUXVb;
    %
    local_FLUXTb = FLUXCb1(iBFace)*pressure(iElement1) + FLUXCb2(iBFace)*pressure_b(k) + FLUXVb(iBFace);
    theFluxes.FLUXTf(iBFaces(iBFace)) = local_FLUXTb;
    mdot_b(iBFace) = mdot_b(iBFace) + local_FLUXTb;
    %
end
theMdotField.phi = mdot_b;
cfdSetMeshField(theMdotField);
end

%===================================================
% INLET-Specified Average Pressure Condition
%===================================================

function theFluxes =  cfdAssembleMdotTermInletSpecifiedAveragePressureBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

theFluxes = cfdAssembleMdotTermInletFixedValueBC(iPatch,theMesh,theCoefficients,theContinuityEquation);
end


%===================================================
% OUTLET-Outlet Condition
%===================================================

function theFluxes =  cfdAssembleMdotTermOutletOutletBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
Sf = -[theMesh.faces(iBFaces).Sf]';
%
% Get Arrays
%
%
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_previous_b = theMdotField.phi(iBFaces);
mdot_b = zeros(size(mdot_previous_b));

theVelocityField = cfdGetMeshField('U');
vel_b = theVelocityField.phi(iBElements,:);

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi(iBElements);

thePressureField = cfdGetMeshField('p');
pressure = thePressureField.phi(iOwners);
pressure_b =thePressureField.phi(iBElements);
%
% ---------
%  STEP 1
% ---------
%  Assemble FLUXCb_1, FLUXCb_2 and FLUXVb coefficients
%
%---------------------------------------------------------------
% assemble RHIE-CHOW Interpolation Term I, II, III, VIII and IX
%---------------------------------------------------------------											%
% loop over patchFaces
fvmFaces = theMesh.faces;
theNumberOfPatchFaces = theMesh.boundaries(iPatch).numberOfBFaces;
patchFacesIndices = find([fvmFaces(:).patchIndex] == iPatch);
%
for k=1:numberOfBFaces
    %
    iBFace = iFaceStart-1+k;
    iOwner = iOwners(k);
    %
    local_FLUXCb1 =  0.0;
    local_FLUXCb2 =  0.0;
    local_FLUXVb =  0.0;
    %
    % ONLY TERM assemble term I
    %
    U_b = Sf(k)*vel_b(k)';
    local_FLUXVb = rho_b(k)*U_b;  % term I
    %
    theFluxes.FLUXC1f(iBFace) = theFluxes.FLUXC1f(iBFace) + local_FLUXCb1;
    theFluxes.FLUXC2f(iBFace) = theFluxes.FLUXC2f(iBFace) + local_FLUXCb2;
    theFluxes.FLUXVf(iBFace)  = theFluxes.FLUXVf(iBFace) + local_FLUXVb;
    %
    % assemble terms X
    %     (mdot_f/rho_f)*(?rho/?p) P'
    %
    %double local_mdot_b  = local_FLUXCb_1*(pressure[iElement1]+ Pref) + local_FLUXCb_2*(pressure_b[iBFace]+ Pref) + local_FLUXVb;
    %
    %local_FLUXCb1 = local_FLUXCb1 + fmax(local_mdot_b/rho_b(iBFace),0.0)*drhodp_b(iBFace);
    %local_FLUXCb2 = local_FLUXCb2 - fmax(-local_mdot_b/rho_b(iBFace),0.0)*drhodp_b(iBFace);
    %local_FLUXVb  = local_FLUXVb -(fmax(local_mdot_b/rho_b(iBFace),0.0)*drhodp_b(iBFace)*(pressure(iElement1)+ Pref) -fmax(-local_mdot_b/rho_b[iBFace],0.0)*drhodp_b[iBFace]*(pressure_b[iBFace]+ Pref));
    % local_FLUXVb += - local_mdot_b;
    %
    local_FLUXTb = theFluxes.FLUXC1f(iBFace)*pressure(k) + theFluxes.FLUXC2f(iBFace)*pressure_b(k) + theFluxes.FLUXVf(iBFace);
    theFluxes.FLUXTf(iBFace) = local_FLUXTb;
    mdot_b(k) = mdot_b(k) + local_FLUXTb;
end
theMdotField.phi(iBFaces) = mdot_b;
cfdSetMeshField(theMdotField);
end



%===================================================
% OUTLET-specifiedValue Condition
%===================================================

function theFluxes =  cfdAssembleMdotTermOutletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iFluid)
%
theMesh = cfdGetMesh;
%
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
%
Sb = [theMesh.faces(iBFaces).Sf]';
CN = [theMesh.faces(iBFaces).CN]';
eCN = [theMesh.faces(iBFaces).eCN]';
%
% Get Arrays
%
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');

theVelocityField = cfdGetMeshField('U');
vel_b = theVelocityField.phi(iBElements,:);

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi(iBElements);

thePressureField = cfdGetMeshField('p');
pressureC = thePressureField.phi(iOwners);
pressure_b = thePressureField.phi(iBElements);
pressure_grad_b = thePressureField.phiGradient(iBElements,:);

%
% Get Arrays
%
theD1Field = cfdGetMeshField('DU1');
theD2Field = cfdGetMeshField('DU2');
theD3Field = cfdGetMeshField('DU3');
%
d1_b = theD1Field.phi(iBElements);
d2_b = theD2Field.phi(iBElements);
d3_b = theD3Field.phi(iBElements);

%
% ---------
%  STEP 1
% ---------
%  Assemble FLUXCb1, FLUXCb2 and FLUXVb coefficients
%
%---------------------------------------------------------------
% assemble RHIE-CHOW Interpolation Term I, II, III, VIII and IX
%---------------------------------------------------------------
%
local_FLUXCb1 =  0.0;
local_FLUXCb2 =  0.0;
local_FLUXVb  =  0.0;
%
% Assemble Coefficients
%
DU1_b = d1_b;
DU2_b = d2_b;
DU3_b = d3_b;

DUSb = [DU1_b.*Sb(:,1),DU2_b.*Sb(:,2),DU3_b.*Sb(:,3)];
eDUSb = [DUSb(:,1)./cfdMagnitude(DUSb),DUSb(:,2)./cfdMagnitude(DUSb),DUSb(:,3)./cfdMagnitude(DUSb)];

DUEb = [cfdMagnitude(DUSb)./dot(eCN(:,:)',eDUSb(:,:)')'.*eCN(:,1),cfdMagnitude(DUSb)./dot(eCN(:,:)',eDUSb(:,:)')'.*eCN(:,2),cfdMagnitude(DUSb)./dot(eCN(:,:)',eDUSb(:,:)')'.*eCN(:,3)];
geoDiff = cfdMagnitude(DUEb)./cfdMagnitude(CN);

DUTb = DUSb - DUEb;
%
% Assemble term I
%
U_bar_b = dot(vel_b(:,:)',Sb(:,:)')';
local_FLUXVb = local_FLUXVb + rho_b.*U_bar_b;  % term I
%
% Assemble term II and linearize it
%
local_FLUXCb1 = local_FLUXCb1 + rho_b.*geoDiff;
local_FLUXCb2 = local_FLUXCb2 - 0.0;
local_FLUXVb  = local_FLUXVb  - rho_b.*dot(pressure_grad_b(:,:)',DUTb(:,:)')'; % Neglect non-orthogonal term
%
% Assemble term III
%
local_FLUXVb  = local_FLUXVb + rho_b.*dot(pressure_grad_b(:,:)',DUSb(:,:)')';
local_FLUXTb = local_FLUXCb1.*pressureC + local_FLUXCb2.*pressure_b + local_FLUXVb;

mdot_b = local_FLUXTb;
%
% Store mdot_b
%
theMdotField.phi(iBFaces) = mdot_b;
cfdSetMeshField(theMdotField);
%
% Assemble into global terms
%
theFluxes.FLUXC1f(iBFaces) = local_FLUXCb1;
theFluxes.FLUXC2f(iBFaces) = local_FLUXCb2;
theFluxes.FLUXVf(iBFaces)  = local_FLUXVb;
theFluxes.FLUXTf(iBFaces) = local_FLUXTb;

end


