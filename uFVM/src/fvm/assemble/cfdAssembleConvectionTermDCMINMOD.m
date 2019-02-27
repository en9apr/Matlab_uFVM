function theCoefficients = cfdAssembleConvectionTermDCMINMOD(theMesh,theCoefficients,theScalarFieldUserName,theTerm)
%===================================================

% This function computed the correction to the UPWIND convection assembly
% to yield the SOU scheme assemble using the deferred correction method
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================

%
% Get scalar name and fluidtag
%
theScalarFieldName = cfdConvertName(theScalarFieldUserName);
theFluidIndex = cfdGetFluidIndex(theScalarFieldName);
theFluidTag = cfdGetFluidTag(theScalarFieldName);
theFieldBaseName = cfdGetBaseName(theScalarFieldName);

%
% Get the RHS
%
b = theCoefficients.RHS;
%
% Get the VF for multiphase case
%
theVFfField = cfdGetModel(['VF_f' theFluidTag]);
vf_f = theVFfField.phi;
%
% Get elements and faces
%
fvmElements = theMesh.elements;
fvmFaces = theMesh.faces;
%
% Get mdot
%
theMdotName = ['Mdot' theFluidTag];
mdotField = cfdGetModel(theMdotName);
mdot_f = mdotField.phi;
%
% Get the face and control volume gradients
%
theScalarField = cfdGetModel(theScalarFieldName);
phiGradients = theScalarField.phiGradient;
phiGradients_f = fvmInterpolateFieldGradientToFaces(theMesh,theScalarField);
%//////////////////////////////////////////////////////////
%
% Loop over all interior faces and assemble the face fllux correction
%
interiorFacesIndices = find([fvmFaces(:).patchIndex]==0);
for iFace=interiorFacesIndices
    theFace = fvmFaces(iFace);
    iElement1 = theFace.element1;
    iElement2 = theFace.element2;
    %
    % Determine the C (Upwind) control volume
    %
    if(mdot_f(iFace)>0)
        iElementC = iElement1;
        iElementD = iElement2;
    else
        iElementC = iElement2;
        iElementD = iElement1;
    end
    % 
    % retrieve the needed values for the C Control Volume
    %
    theElementC = fvmElements(iElementC);
    theElementD = fvmElements(iElementD);
    r_CD = theElementD.centroid - theElementC.centroid;
    phiGradientC = phiGradients(iElementC,:);
    phiC= theScalarField.phi(iElementC);
    phiD= theScalarField.phi(iElementD);
    phiMin = min(theScalarField.phi([theElementC.neighbours iElementC]));
    phiMax = max(theScalarField.phi([theElementC.neighbours iElementC]));
    phiU = phiD -2*phiGradientC*r_CD';
    phiU=max(min(phiMax,phiU),phiMin);
    %
    % Compute phiNormC
    %
       %
    % compute phiNormf
    %   
    if(abs(phiD-phiU)>1e-9)
        phiNormC = (phiC-phiU)/(phiD-phiU);
        phiNormf = phiNormC;
        if((phiNormC>=0) && (phiNormC<=1)) 
             phiNormf = min(1.5*phiNormC,0.5*phiNormC+0.5);
        end
        phiHRf = phiNormf*(phiD-phiU)+phiU;         
    else
        phiHRf = phiC;
    end
    %
    % Compute the correction flux
    % -mdot*(phi_f_HR - phi_f_UPWIND)
    %
    FLUXTf = -mdot_f(iFace)*(phiHRf - phiC);
    %
    % assemble into RHS
    %
    b(iElement1) = b(iElement1) + FLUXTf*vf_f(iFace);
    
    b(iElement2) = b(iElement2) - FLUXTf*vf_f(iFace);
end 


theCoefficients.RHS = b;
