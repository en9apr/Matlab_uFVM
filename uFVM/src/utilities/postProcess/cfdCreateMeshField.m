function theMeshField = cfdCreateMeshField(theFieldName, theFoamField)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function creates a mesh field from FOAM data
%--------------------------------------------------------------------------

% Get mesh info
theMesh = cfdGetMesh;
theFaces = theMesh.faces;
theElements = theMesh.elements;

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

% Define mesh field
if strcmp(theFoamField.class,'volScalarField')
    locale = 'Elements';
elseif strcmp(theFoamField.class,'volVectorField')
    locale = 'Elements';
elseif strcmp(theFoamField.class,'surfaceScalarField')
    locale = 'Faces';
end

% Define mesh field
theMeshField = cfdGetMeshField(theFieldName, locale);

% Evaluate internal field
internalField = theFoamField.internalField;
value = internalField.value;
if strcmp(internalField.valueType, 'uniform')
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
        theMeshField.phi(1:numberOfElements, 1) = value;
    elseif strcmp(theFoamField.class, 'volVectorField')
        for iComponent=1:3
            theMeshField.phi(1:numberOfElements, iComponent) = value(:, iComponent);
        end
    elseif strcmp(theFoamField.class, 'surfaceScalarField')
        theMeshField.phi(1:numberOfInteriorFaces, 1) = value;
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
cfdSetMeshField(theMeshField);

% Create model except when the field is the mdot_f field
if ~strcmp(theFieldName, 'mdot_f')
    theFieldModel.userName = theFieldName;
    theFieldModel.name = theFieldName;
    theFieldModel.class = 'Constant';
    theFieldModel.type = theMeshField.type;
    theFieldModel.urf = 1.0;
    theFieldModel.calcGradient = false;
    theFieldModel.gradientType = '';
    
    if theFieldModel.calcGradient
        theFieldModel.phiGradient = [];
    end
    
    cfdSetModel(theFieldModel);
end

end

