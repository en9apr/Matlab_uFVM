function cfdUpdateTermCoefficientField(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function updates the term coefficient fields
%--------------------------------------------------------------------------

theField = cfdGetModel(theEquationName);
theTerms = theField.terms;

theMesh = cfdGetMesh;

if strcmp(theField.name, 'U')
    for iTerm=1:length(theTerms)
        fmvOperator = theTerms{iTerm}.fvmOperator;
        if strcmp(fmvOperator, 'ddt')
            theRhoField = cfdGetMeshField(['rho_',theField.name,'eq']);
            theRhoExpression = theField.rhoName;
            
            if ~isempty(theRhoExpression)
                % Split the rho expression into its constituting variables
                coefficientVariables = retrieveVariables(theRhoExpression);
                for iVariable=1:length(coefficientVariables)
                    eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                    eval([coefficientVariables{iVariable},' = theVariableField.phi;']);
                end
                
                % Adjust expression to include vectorial arithmetic operators
                theRhoExpression = strrep(theRhoExpression, '*', '.*');
                theRhoExpression = strrep(theRhoExpression, '/', './');
                
                % Update the rho field of U equation
                theRhoField.phi = eval(theRhoExpression);
            else
                % Set unity rho field
                theRhoField.phi = ones(size(theRhoField.phi));
            end
            
            % Store in data base
            cfdSetMeshField(theRhoField);
        elseif strcmp(fmvOperator, 'laplacian')
            theGammaField = cfdGetMeshField(['gamma_',theField.name,'eq'], 'Faces');
            theGammaExpression = theField.gammaName;
            
            if ~isempty(theGammaExpression)
                % Split the gamma expression into its constituting variables.
                % And interpolate variables to faces
                coefficientVariables = identifyFields(theGammaExpression);
                for iVariable=1:length(coefficientVariables)
                    eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                    eval('phi_f = cfdInterpolateFromElementsToFaces(''Average'', theVariableField.phi);');
                    eval([coefficientVariables{iVariable},' = phi_f;']);
                end
                
                % Adjust expression to include vectorial arithmetic operators
                theGammaExpression = strrep(theGammaExpression, '*', '.*');
                theGammaExpression = strrep(theGammaExpression, '/', './');
                
                % Update the gamma field of U equation
                theGammaField.phi = eval(theGammaExpression);
            else
                % Set unity gamma field
                theGammaField.phi = ones(size(theGammaField.phi));
            end
            
            % Store in data base
            cfdSetMeshField(theGammaField);
        end
    end
    
    % Update the rho field if not yet (this usually happens when transient
    % term is not included in the U equation)
    theTerm = cfdGetTermInEquation('U', 'Transient');
    if isempty(theTerm)
        theRhoField = cfdGetMeshField('rho_Ueq');
        
        % The following is only applicable if a convection term exists
        if isfield(theField, 'psiName')
            thePsiExpression = theField.psiName;
            
            % Get the rho expression but removing the vector U
            theRhoExpression = strrep(thePsiExpression, 'U', '1');
            
            coefficientVariables = retrieveVariables(theRhoExpression);
            for iVariable=1:length(coefficientVariables)
                eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                eval([coefficientVariables{iVariable},' = theVariableField.phi;']);
            end
            
            % Adjust expression to include vectorial arithmetic operators
            theRhoExpression = strrep(theRhoExpression, '*', '.*');
            theRhoExpression = strrep(theRhoExpression, '/', './');
            
            % Update the psi field of U equation
            if strcmp(theRhoExpression, '1')
                theRhoField.phi = ones(size(theRhoField.phi));
            else
                theRhoField.phi = eval(theRhoExpression);
            end
            
            cfdSetMeshField(theRhoField);
        else
            theRhoField = cfdGetMeshField('rho');
            if isempty(theRhoField)
                error('\n%s\n', 'You have to define transport property rho to proceed');
            end
        end
    end
    
elseif ~strcmp(theField.name, 'p')
    for iTerm=1:length(theTerms)
        fmvOperator = theTerms{iTerm}.fvmOperator;
        if strcmp(fmvOperator, 'ddt')
            theRhoField = cfdGetMeshField(['rho_',theField.name,'eq']);
            theRhoExpression = theField.rhoName;
            
            if ~isempty(theRhoExpression)
                % Split the rho expression into its constituting variables
                coefficientVariables = retrieveVariables(theRhoExpression);
                for iVariable=1:length(coefficientVariables)
                    eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                    eval([coefficientVariables{iVariable},' = theVariableField.phi;']);
                end
                
                % Adjust expression to include vectorial arithmetic operators
                theRhoExpression = strrep(theRhoExpression, '*', '.*');
                theRhoExpression = strrep(theRhoExpression, '/', './');
                
                % Update the rho field of U equation
                theRhoField.phi = eval(theRhoExpression);
            else
                % Set unity rho field
                theRhoField.phi = ones(size(theRhoField.phi));
            end
            
            % Store in data base
            cfdSetMeshField(theRhoField);
        elseif strcmp(fmvOperator, 'div')
            % Check if U field exists because this is a convection term
            % which is usually associated with U field
            theVelocityField = cfdGetMeshField('U');
            if isempty(theVelocityField)
                error('\n%s\n', 'U field is not defined. You have to define a U equation, or you may assign U as a property.');
            end
            
            thePsiField = cfdGetMeshField(['psi_',theField.name,'eq'], 'Faces');
            thePsiExpression = theField.psiName;
            
            % Check if U is included in the psi field
            if ~strfind(thePsiExpression, 'U')
                error('\n%s\n', ['Check the ',theField.name,' equation. In the div(psi, ',theField.name,'), psi must include U.']);
            end
            
            % Split the psi expression into its constituting variables.
            % and interpolate variables to faces
            theModPsiExpression = strrep(thePsiExpression, 'U', '1');
            coefficientVariables = retrieveVariables(theModPsiExpression);
            for iVariable=1:length(coefficientVariables)
                eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                eval('phi_f = cfdInterpolateFromElementsToFaces(''Average'', theVariableField.phi);');
                eval([coefficientVariables{iVariable},' = phi_f;']);
            end
            
            % Adjust expression to include vectorial arithmetic operators
            theModPsiExpression = strrep(theModPsiExpression, '*', '.*');
            theModPsiExpression = strrep(theModPsiExpression, '/', './');
            coefficient = eval(theModPsiExpression);
            
            U = theVelocityField.phi;
            U_f = cfdInterpolateFromElementsToFaces('Average', U);
            
            iFaces = 1:theMesh.numberOfFaces;
            Sf = [theMesh.faces(iFaces).Sf]';
            DeltaVol = [theMesh.faces(iFaces).DeltaVol]';
            
            if cfdIsTransient
                dt = cfdGetDt;
                
                % Initialize psi at interior faces
                psi1 = coefficient .* dot(Sf',U_f')';
                psi2 = coefficient .* DeltaVol / dt;
                thePsiField.phi = psi1 - psi2;
            else
                % Update the psi field of U equation
                iFaces = 1:theMesh.numberOfFaces;
                Sf = [theMesh.faces(iFaces).Sf]';
                thePsiField.phi = coefficient .* dot(Sf',U_f')';
            end
            
            % Store in data base
            cfdSetMeshField(thePsiField);
        elseif strcmp(fmvOperator, 'laplacian')
            theGammaField = cfdGetMeshField(['gamma_',theField.name,'eq'], 'Faces');
            theGammaExpression = theField.gammaName;
            
            if ~isempty(theGammaExpression)
                % Split the gamma expression into its constituting variables.
                % And interpolate variables to faces
                coefficientVariables = retrieveVariables(theGammaExpression);
                for iVariable=1:length(coefficientVariables)
                    eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                    eval('phi_f = cfdInterpolateFromElementsToFaces(''Average'', theVariableField.phi);');
                    eval([coefficientVariables{iVariable},' = phi_f;']);
                end
                
                % Adjust expression to include vectorial arithmetic operators
                theGammaExpression = strrep(theGammaExpression, '*', '.*');
                theGammaExpression = strrep(theGammaExpression, '/', './');
                
                % Update the gamma field of U equation
                theGammaField.phi = eval(theGammaExpression);
                
            else
                % Set unity gamma field
                theGammaField.phi = ones(size(theGammaField.phi));
            end
            
            % Store in data base
            cfdSetMeshField(theGammaField);
        end
    end
    
    % Update the rho field if not yet (this usually happens when transient
    % term is not included in the U equation)
    theTerm = cfdGetTermInEquation(theField.name, 'Transient');
    if isempty(theTerm)
        theRhoField = cfdGetMeshField(['rho_',theField.name,'eq']);
        
        % The following is only applicable if a convection term exists
        if isfield(theField, 'psiName')
            thePsiExpression = theField.psiName;
            
            % Get the rho expression but removing the vector U
            theRhoExpression = strrep(thePsiExpression, 'U', '1');
            
            coefficientVariables = retrieveVariables(theRhoExpression);
            for iVariable=1:length(coefficientVariables)
                eval(['theVariableField = cfdGetMeshField(''',coefficientVariables{iVariable},''');']);
                eval([coefficientVariables{iVariable},' = theVariableField.phi;']);
            end
            
            % Adjust expression to include vectorial arithmetic operators
            theRhoExpression = strrep(theRhoExpression, '*', '.*');
            theRhoExpression = strrep(theRhoExpression, '/', './');
            
            % Update the psi field of U equation
            if strcmp(theRhoExpression, '1')
                theRhoField.phi = ones(size(theRhoField.phi));
            else
                theRhoField.phi = eval(theRhoExpression);
            end
            
            cfdSetMeshField(theRhoField);
        else
%             theRhoField = cfdGetMeshField('rho');
%             if isempty(theRhoField)
%                 error('\n%s\n', 'You have to define transport property rho to proceed');
%             end
        end
    end
end
