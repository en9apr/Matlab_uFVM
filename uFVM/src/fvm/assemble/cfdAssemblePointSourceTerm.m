function theCoefficients = fvmAssemblePointSourceTerm(theCoefficients,theScalarFieldUserName,theTerm)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
theMesh = cfdGetMesh;

theScalarField = cfdGetModel(theScalarFieldUserName);


myPoint = theTerm.position;
Sp = theTerm.Sp;

%
theNumberOfElements = theMesh.numberOfElements;
%
b = theCoefficients.RHS;

%
% Assemble Source Term
%
for iElement=1:theNumberOfElements
    %
    % Compute Fluxes
    %
    if(isPointInElement(myPoint,iElement))
         b(iElement)  = b(iElement) + Sp;
    end
end
%
theCoefficients.RHS = b;

