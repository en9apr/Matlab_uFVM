function theMeshField = createMeshField(theFieldName, theFoamField)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function creates a mesh field from FOAM data
%--------------------------------------------------------------------------

% Define mesh field
theMeshField.userName = theFieldName;
theMeshField.name = theFieldName;
if strcmp(theFoamField.class,'volScalarField')
    theMeshField.type = 'Scalar';
    theMeshField.locale = 'Elements';
elseif strcmp(theFoamField.class,'volVectorField')
    theMeshField.type = 'Vector';
    theMeshField.locale = 'Elements';
elseif strcmp(theFoamField.class,'surfaceScalarField')
    theMeshField.type = 'Scalar';
    theMeshField.locale = 'Faces';
end

theMesh = cfdGetMesh;

numberOfElements = theMesh.numberOfElements;
numberOfBElements = theMesh.numberOfBElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theMesh.numberOfBFaces;

% Initialize phi and phiGradient
if strcmp(theFoamField.class,'volScalarField')
    theMeshField.phi = zeros(numberOfElements + numberOfBElements, 1);
    theMeshField.phiGradient = zeros(numberOfElements + numberOfBElements, 3, 1);
elseif strcmp(theFoamField.class,'volVectorField')
    theMeshField.phi = zeros(numberOfElements + numberOfBElements, 3);
    theMeshField.phiGradient = zeros(numberOfElements + numberOfBElements, 3, 3);
elseif strcmp(theFoamField.class,'surfaceScalarField')
    theMeshField.phi = zeros(numberOfInteriorFaces + numberOfBFaces, 1);
end

% Evaluate internal field
internalField = theFoamField.internalField;
if strcmp(internalField.valueType, 'uniform')
    value = internalField.value;
    if strcmp(theFoamField.class, 'volScalarField')
        theMeshField.phi(1:numberOfElements, 1) = value .* ones(numberOfElements, 1);
    elseif strcmp(theFoamField.class, 'volVectorField')
        for iComponent=1:3
            theMeshField.phi(1:numberOfElements, iComponent) = value(iComponent) .* ones(numberOfElements, 1);
        end
    elseif strcmp(theFoamField.class, 'surfaceScalarField')
        theMeshField.phi(1:numberOfInteriorFaces, 1) = value .* ones(numberOfInteriorFaces, 1);
    end
elseif strcmp(internalField.valueType, 'nonuniform')
    if strcmp(theFoamField.class, 'volScalarField')
        theMeshField.phi(1:numberOfElements, 1) = internalField.value;
    elseif strcmp(theFoamField.class, 'volVectorField')
        for iComponent=1:3
            theMeshField.phi(1:numberOfElements, iComponent) = internalField.value(:, iComponent);
        end
    elseif strcmp(theFoamField.class, 'surfaceScalarField')
        theMeshField.phi(1:numberOfInteriorFaces, 1) = internalField.value;       
    end        
end

% Evaluate boundary field
theNumberOfPatches = theMesh.numberOfBoundaries;
phi = theMeshField.phi;
for iPatch=1:theNumberOfPatches           
    
    % Boundary Mesh info
    theBoundary = theMesh.boundaries(iPatch);
    %
    numberOfBFaces = theBoundary.numberOfBFaces;
    iFaceStart = theBoundary.startFace;
    iFaceEnd = iFaceStart+numberOfBFaces-1;
    iBFaces = iFaceStart:iFaceEnd;
    %
    iBElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
    iBElementEnd = iBElementStart+numberOfBFaces-1;
    iBElements = iBElementStart:iBElementEnd;
    %    
    iOwners = [theMesh.faces(iBFaces).iOwner];
    %
    theLocale = ['BPatch' num2str(iPatch)];
    
    % Get boundary field
    boundaryField = theFoamField.boundaryField{iPatch};    
    
    if strcmp(boundaryField.type, 'fixedValue') || strcmp(boundaryField.type, 'calculated')       
        if strcmp(boundaryField.valueType, 'uniform')
            if strcmp(theFoamField.class, 'volVectorField')
                value = ['[',num2str(boundaryField.value(1)),';',num2str(boundaryField.value(2)),';',num2str(boundaryField.value(3)),']'];
                phi_b = cfdComputeFormulaAtLocale(value,theLocale,'Vector');
                phi(iBElements, :) = phi_b;                
            elseif strcmp(theFoamField.class, 'volScalarField')
                value = num2str(boundaryField.value);
                phi_b = cfdComputeFormulaAtLocale(value,theLocale,'Scalar');
                phi(iBElements) = phi_b;
            elseif strcmp(theFoamField.class, 'surfaceScalarField')
                value = num2str(boundaryField.value); 
                phi_b = cfdComputeFormulaAtLocale(value,theLocale,'Scalar');
                phi(iBFaces) = phi_b;
            end            
        elseif strcmp(boundaryField.valueType, 'nonuniform')
            value = boundaryField.value;
            if strcmp(theFoamField.class, 'volVectorField')                
                phi(iBElements, :) = value;                
            elseif strcmp(theFoamField.class, 'volScalarField')
                phi(iBElements) = value;
            elseif strcmp(theFoamField.class, 'surfaceScalarField')
                phi(iBFaces) = value;
            end             
        end        
    elseif strcmp(boundaryField.type, 'noSlip')
        if strcmp(theFoamField.class, 'volVectorField')
            value = '[0;0;0]';
            phi_b = cfdComputeFormulaAtLocale(value,theLocale,'Vector');
            phi(iBElements, :) = phi_b;
        elseif strcmp(theFoamField.class, 'volScalarField')
            value = '0';
            phi_b = cfdComputeFormulaAtLocale(value,theLocale,'Scalar');
            phi(iBElements) = phi_b;
        end 
        %
        %
        %        
    elseif strcmp(boundaryField.type, 'fixedGradient')
        d = [theFaces(iBFaces).centroid]' - [theElements(iOwners).centroid]';
        dmag = cfdMagnitude(d);
        e = d./dmag;
        
        % Compute grad phi
        normGrad_b = cfdComputeFormulaAtLocale(value, theLocale);
        n = [theFaces(iBFaces).Sf]' ./ [theFaces(iBFaces).area]';
        grad_b = normGrad_b .* n;
        
        % Evaluate phi at boundary faces
        phi(iBElements) = dmag .* dot(grad_b', e')' + phi(iOwners);
        %
        %
        %
    elseif strcmp(boundaryField.type, 'zeroGradient')
        if strcmp(theFoamField.class, 'volVectorField')
            for iComponent=1:3
                phi(iBElements,iComponent) = phi(iOwners,iComponent);
            end
        elseif strcmp(theFoamField.class, 'volScalarField')
            phi(iBElements) = phi(iOwners);
        end      
        %
        %
        %         
    elseif strcmp(boundaryField.type, 'empty') || strcmp(boundaryField.type, 'symmetry') || strcmp(boundaryField.type, 'slip')
        if strcmp(theFoamField.class, 'volVectorField')
            Sb = [theMesh.faces(iBFaces).Sf]';
            normSb = cfdMagnitude(Sb);
            n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];
            
            phi_normal = dot(phi(iOwners,:)',n')';
            
            for iComponent=1:3
                phi(iBElements,iComponent) = phi(iOwners,iComponent) - phi_normal .* n(:,iComponent);
            end
        elseif strcmp(theFoamField.class, 'volScalarField')
            phi(iBElements) = phi(iOwners);
        end      
        %
        %
        %          
    end             
end
theMeshField.phi = phi;

% Evaluate gradient based on green gauss method
if ~strcmp(theFoamField.class, 'surfaceScalarField')
    theMeshField.phiGradient = cfdComputeGradientGauss0(theMeshField.phi, theMesh);
end