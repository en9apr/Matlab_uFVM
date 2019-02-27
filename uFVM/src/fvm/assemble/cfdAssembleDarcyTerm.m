function theCoefficients = cfdAssembleDarcyTerm(theCoefficients,theScalarFieldUserName,theTerm)
%===================================================
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
theMesh = cfdGetMesh;
fvmElements = theMesh.elements;

theVelocityField = cfdGetModel(theScalarFieldUserName);
theFluidIndex = cfdGetFluidIndex(theScalarFieldUserName);
theFluidTag = cfdGetFluidTag(theScalarFieldUserName);
theFieldBaseName = cfdGetBaseName(theScalarFieldUserName);


%
a = theCoefficients.LHS;
b = theCoefficients.RHS;
%
% x- or y- velocity equation
if(strcmp(theFieldBaseName,'Velx')) 
    e = [1;0];
elseif(strcmp(theFieldBaseName,'Vely')) 
    e = [0;1];
end
% VF
% theVolumeFractionField = cfdGetModel('VF_LIQUID');
% vf = theVolumeFractionField.phi;
%
velc = theVelocityField.phi;
% Kref
Kref = cfdGetConstant('DarcyConstant');

%

theVFField = cfdGetModel(['VF' theFluidTag]);
vf = theVFField.phi;

theNumberOfElements = theMesh.numberOfElements;
for iElement = 1:theNumberOfElements
    volume = fvmElements(iElement).volume;  
    %
    % Compute Fluxes
    %
    FLUXCE =  volume * Kref *(exp(1-vf(iElement))-1);
    FLUXVE = 0;
    FLUXTE = FLUXCE*velc(iElement) + FLUXVE;
    %
    % Assemble Fluxes
    %
    a(iElement,iElement) =  a(iElement,iElement) + FLUXCE*vf(iElement);    
    b(iElement) = b(iElement) - FLUXTE*vf(iElement);
    
end
%
theCoefficients.LHS = a;
theCoefficients.RHS = b;
