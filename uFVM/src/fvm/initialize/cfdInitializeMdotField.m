function cfdInitializeMdotField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes mdot_f term
%--------------------------------------------------------------------------

% get mesh and field info
theMesh = cfdGetMesh;
theEquation = cfdGetModel('U');
theVelocityField = cfdGetMeshField('U');
theMdotField = cfdGetMeshField('mdot_f', 'Faces');

% Retrieve the density expression from the psi field
thePsiExpression = theEquation.psiName;

if ~isempty(thePsiExpression)
    %
    % If divergence term exists
    %
    
    % Extract density expression
    theDensityExpression = strrep(thePsiExpression, 'U', '1');
    
    if ~strcmp(theDensityExpression, '1')
        % Check if U variable exist in the psi expression. If nothing change in the
        % above command, then U hasn't been there!
        if strcmp(theDensityExpression, thePsiExpression)
            error('\n%s\n', 'Check the U equation. In the div(psi, U), psi must include vector U.');
        end
        
        % Split the density expression into its constituting variables
        coefficientVariables = retrieveVariables(theDensityExpression);
        for iVariable=1:length(coefficientVariables)
            eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
            eval([coefficientVariables{iVariable},' = theVariableField.phi;']);
        end
        
        % Adjust expression to include vectorial arithmetic operators
        theDensityExpression = strrep(theDensityExpression, '*', '.*');
        theDensityExpression = strrep(theDensityExpression, '/', './');
        
        % Update the density field of U equation
        rho = eval(theDensityExpression);
    else
        % Set unity gamma field
        rho = ones(size(theMdotField.phi));
    end
    
    if cfdIsTransient
        % Get required fields
        dt = cfdGetDt;
        U = theVelocityField.phi;
        
        % interpolate fields to faces
        U_f = cfdInterpolateFromElementsToFaces('Average', U);
        rho_f = cfdInterpolateFromElementsToFaces('Average', rho);
        
        % Initialize mdot_f at interior faces
        iFaces = 1:theMesh.numberOfFaces;
        Sf = [theMesh.faces(iFaces).Sf]';
        DeltaVol = [theMesh.faces(iFaces).DeltaVol]';
        mdot_f1 = rho_f .* dot(Sf',U_f')';
        mdot_f2 = rho_f .* DeltaVol / dt;
        mdot_f = mdot_f1 - mdot_f2;
        
        theMdotField.phi = mdot_f;
        cfdSetMeshField(theMdotField);
    else
        U = theVelocityField.phi;
        
        % interpolate fields to faces
        U_f = cfdInterpolateFromElementsToFaces('Average',U);
        rho_f = cfdInterpolateFromElementsToFaces('Average',rho);
        
        % Initialize mdot_f at interior faces
        iFaces = 1:theMesh.numberOfFaces;
        Sf = [theMesh.faces(iFaces).Sf]';
        mdot_f = rho_f.*dot(Sf',U_f')';
        
        theMdotField.phi = mdot_f;
        cfdSetMeshField(theMdotField);
    end
    
else
    %
    % If divergence term doesn't exist
    %
    theDensityField = cfdGetMeshField('rho');
    if cfdIsTransient
        % Get required fields
        dt = cfdGetDt;
        U = theVelocityField.phi;
        rho = theDensityField.phi;
        
        % interpolate fields to faces
        U_f = cfdInterpolateFromElementsToFaces('Average', U);
        rho_f = cfdInterpolateFromElementsToFaces('Average', rho);
        
        % Initialize mdot_f at interior faces
        iFaces = 1:theMesh.numberOfFaces;
        Sf = [theMesh.faces(iFaces).Sf]';
        DeltaVol = [theMesh.faces(iFaces).DeltaVol]';
        mdot_f1 = rho_f .* dot(Sf',U_f')';
        mdot_f2 = rho_f .* DeltaVol / dt;
        mdot_f = mdot_f1 - mdot_f2;
        
        theMdotField.phi = mdot_f;
        cfdSetMeshField(theMdotField);
    else
        U = theVelocityField.phi;
        rho = theDensityField.phi;
        
        % interpolate fields to faces
        U_f = cfdInterpolateFromElementsToFaces('Average',U);
        rho_f = cfdInterpolateFromElementsToFaces('Average',rho);
        
        % Initialize mdot_f at interior faces
        iFaces = 1:theMesh.numberOfFaces;
        Sf = [theMesh.faces(iFaces).Sf]';
        mdot_f = rho_f.*dot(Sf',U_f')';
        
        theMdotField.phi = mdot_f;
        cfdSetMeshField(theMdotField);
    end
    
end

% Initialize mdot field for each fluid if the application class is
% multiphase
applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'multiphase')
    theNumberOfFluids = cfdGetNumberOfFluids;
    for iFluid=1:theNumberOfFluids
        theMdotFieldName = ['mdot_f', num2str(iFluid)];
        theMdotField = cfdGetMeshField(theMdotFieldName, 'Faces');
        theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
        if cfdIsTransient
            % Get required fields
            dt = cfdGetDt;
            U = theVelocityField.phi;
            rho = theDensityField.phi;
            
            % Interpolate fields to faces
            U_f = cfdInterpolateFromElementsToFaces('Average', U);
            rho_f = cfdInterpolateFromElementsToFaces('Average', rho);
            
            % Initialize mdot_f at interior faces
            iFaces = 1:theMesh.numberOfFaces;
            Sf = [theMesh.faces(iFaces).Sf]';
            DeltaVol = [theMesh.faces(iFaces).DeltaVol]';
            mdot_f1 = rho_f .* dot(Sf',U_f')';
            mdot_f2 = rho_f .* DeltaVol / dt;
            mdot_f = mdot_f1 - mdot_f2;
            
            theMdotField.phi = mdot_f;
            cfdSetMeshField(theMdotField);
        else
            U = theVelocityField.phi;
            rho = theDensityField.phi;
            
            % interpolate fields to faces
            U_f = cfdInterpolateFromElementsToFaces('Average',U);
            rho_f = cfdInterpolateFromElementsToFaces('Average',rho);
            
            % Initialize mdot_f at interior faces
            iFaces = 1:theMesh.numberOfFaces;
            Sf = [theMesh.faces(iFaces).Sf]';
            mdot_f = rho_f.*dot(Sf',U_f')';
            
            theMdotField.phi = mdot_f;
            cfdSetMeshField(theMdotField);
        end
    end
end
