function cfdDefineEquation(theEquationName, theEquationExpression)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% cfdSetupEquation   Setup cfd equation.
%     cfdSetupEquation(theEquationName,theEquationExpression) creates an
%     equation of name theEquationName. Additional and often necessary
%     argument can be added for the expression of the equation.
%
%     Examples:
%
%     % Create incompressible momentum equation U with convection,
%     % diffusion and pressure gradient terms:
%     cfdSetupEquation('U','div(rho*U*U)=laplacian(mu*U)-grad(p)')
%
%     % Create pressure equation p:
%     cfdSetupEquation('p')
%
%     % Create incompressible temperature equation T:
%     theEquation = cfdSetupEquation('T','div(rho*Cp*U*U)=laplacian(k*T)')
%--------------------------------------------------------------------------

% Application class
applicationClass = cfdGetApplicationClass;

if nargin==1 && strcmp(theEquationName, 'p')
    if strcmp(applicationClass, 'incompressible')
        theEquationExpression = 'div(U) = 0';
    elseif strcmp(applicationClass, 'compressible')
        theEquationExpression = 'ddt(rho) + div(rho*U) = 0';
    end    
end

% Convert equation name in case if availablilty of phases. No changes are
% made to the name in case of single phase flow
theEquationUserName = cfdConvertFieldPhaseName(theEquationName);

% Check if a file of the same name of the equation name exists
if exist(['0/',theEquationUserName], 'file')~=2
    list = dir;
    folderName = list(1).folder;
    [pathstr,name,ext] = fileparts(folderName);
    error('\n%s\n', ['A file of the name ', theEquationUserName, ' doesn''t exist in ', name,'/0 directory.']);
end

% Print Equation
fprintf('\n-------------------------------------------------');
fprintf(['\nSolving ', theEquationName,' Eq:\n', theEquationExpression,'\n']);

% Call data base
global Domain;

% Create theEquation structure and store its name and username
theEquation.name = theEquationName;
theEquation.userName = theEquationUserName;
theEquation.class = 'Equation';

% Retrieve the fields defined by foam
theField = getfield(Domain.foam.fields, theEquationName);

% The Equation Type
theFieldType = theField.class;
if strcmp(theFieldType, 'volScalarField')
    theEquation.type = 'Scalar';
elseif strcmp(theFieldType, 'volVectorField')
    theEquation.type = 'Vector';
end

% The Equation Initial Conditions
internalField = theField.internalField;
if strcmp(internalField.valueType, 'uniform')
    if strcmp(theFieldType, 'volVectorField')
        theEquation.ic = ['[',num2str(internalField.value(1)),';',num2str(internalField.value(2)), ...
            ';',num2str(internalField.value(3)),']'];
    elseif strcmp(theFieldType, 'volScalarField')
        theEquation.ic = num2str(internalField.value);
    end
elseif strcmp(internalField.valueType, 'nonuniform')
    theEquation.ic = theField.internalField.value;
end

% Under Relaxation Factors
if ~isempty(Domain.foam.fvSolution.relaxationFactors.equations)
    urfFields = Domain.foam.fvSolution.relaxationFactors.fields;
    urfEquations = Domain.foam.fvSolution.relaxationFactors.equations;
    for iField=1:size(urfFields,1)
        if strcmp(urfFields{iField, 1}, theEquationUserName)
            theEquation.urf = str2double(urfFields{iField, 2});
        end
    end
    for iEquation=1:size(urfEquations,1)
        if strcmp(urfEquations{iEquation, 1}, theEquationUserName)
            theEquation.urf = str2double(urfEquations{iEquation, 2});
        end
    end
else
    theEquation.urf = 0.7;
end

% Get the grad scheme of the field
gradSchemes = Domain.foam.fvSchemes.gradSchemes;
theEquation.gradientType = 'Gauss linear';
for iGradScheme=1:size(gradSchemes, 1)
    if strcmp(strrep(['grad(',theEquationName,')'], ' ', ''), ...
            strrep(gradSchemes{iGradScheme, 1}, ' ', ''))
        theEquation.gradientType = gradSchemes{iGradScheme, 2};
        break;
    end
end

% Set other definitions
theEquation.terms = {};
theEquation.residuals = [];
theEquation.source = '';

% Define default AMG options for pressure equation. While set direct
% iterative smoother for other equations
if strcmp(theEquationName,'p')
    theMultigridSolver.isActive = true;
    theMultigridSolver.cycleType = 'W-Cycle';
    theMultigridSolver.maxCycles = '30';
    theMultigridSolver.termination = '0.1';
    theMultigridSolver.maxCoarseLevels = '10';
    theMultigridSolver.preSweep = '1';
    theMultigridSolver.postSweep = '3';
    theMultigridSolver.nCycles = '30';
else
    theMultigridSolver.isActive = false;
end
theEquation.multigrid = theMultigridSolver;

% Set default smoother for all equations
theSolutionMethods.smootherType = 'ILU';
theEquation.theSolutionMethods = theSolutionMethods;

% Save equation/model settings in database
cfdSetupMeshField(theEquationName,'Elements',theEquation.type);
cfdSetModel(theEquation);

% Setup Equation Terms
if strcmp(theEquationName, 'p')
    % Read the equation terms
    %     theEquationTerms = cfdReadEquationTerms(theEquationName, ...
    %         theEquationExpression);
    
    if strcmp(applicationClass, 'multiphase')
        theEquation = cfdAddTerm(theEquationName, 'mdot_mixture_f');
    else
        theEquation = cfdAddTerm(theEquationName, 'mdot_f');
    end         
    
    % Setup DU and DUT fields required for pressure equation assembly
    cfdSetupDU;
    cfdSetupDUT;
else
    % Read equation terms
    theEquationTerms = cfdReadEquationTerms(theEquationName, ...
        theEquationExpression);
    %
    % U Equation
    %
    if strcmp(theEquationName, 'U')
        % Store the details of each term in the equation data base
        theNumberOfTerms = length(theEquationTerms);
        for iTerm=1:theNumberOfTerms
            theTerm = theEquationTerms{iTerm};
            if ~strcmp(theTerm.name, '0') && ~strcmp(theTerm.name, 'Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'fvmOperator', theTerm.operator, ...
                    'coefficientName', theTerm.coefficientName, ...
                    'variableName', theTerm.variableName, ...
                    'scheme', theTerm.scheme, ...
                    'sign', theTerm.sign);
                
                % Setup coefficient fields (excluding divergence term because
                % in this case, the coefficient of the divergence term is
                % mdot_f, which is anyways defined below
                if ~strcmp(theTerm.operator, 'div')
                    cfdSetupTermCoefficientField(theEquationName, theTerm.operator);
                end
            elseif strcmp(theTerm.name, 'Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'sign', theTerm.sign, ...
                    'formula', theTerm.expression, ...
                    'operators', theTerm.operators);
            elseif strcmp(theTerm.name, 'Implicit Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'sign', theTerm.sign, ...
                    'formula', theTerm.expression);               
            end
            cfdPrintTermAttributes(theTerm);
        end
        
        % Check if transient term is included. If not, the rho field
        % (usually associated with the transient term) has to be created on
        % its own
        theTerm = cfdGetTermInEquation('U', 'Transient');
        if isempty(theTerm)
            cfdSetupMeshField('rho_Ueq');
            
            % Update time
            time = cfdGetTime;
            time.type = 'Steady';
            cfdSetTime(time);
        end
        
        % Setup mdot_f field
        theScalarFieldName = 'mdot_f';
        cfdSetupMeshField(theScalarFieldName, 'Faces');
        
        % Setup mdot_f for each fluid if the application class is
        % multiphase
        if strcmp(applicationClass, 'multiphase')
            theNumberOfFluids = cfdGetNumberOfFluids;
            for iFluid=1:theNumberOfFluids
                theScalarFieldName = ['mdot_f', num2str(iFluid)];
                cfdSetupMeshField(theScalarFieldName, 'Faces');
            end
        end
        
        %
        % General (phi) Equation
        %
    else
        % Store the details of each term in the equation data base
        theNumberOfTerms = length(theEquationTerms);
        for iTerm=1:theNumberOfTerms
            theTerm = theEquationTerms{iTerm};
            if ~strcmp(theTerm.name, '0') && ~strcmp(theTerm.name, 'Source') && ~strcmp(theTerm.name, 'Implicit Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'fvmOperator', theTerm.operator, ...
                    'coefficientName', theTerm.coefficientName, ...
                    'variableName', theTerm.variableName, ...
                    'scheme', theTerm.scheme, ...
                    'sign', theTerm.sign);
                
                % Setup coefficient fields
                if ~isempty(theEquationTerms{iTerm}.coefficientName)
                    cfdSetupTermCoefficientField(theEquationName, theTerm.operator);
                end
            elseif strcmp(theTerm.name, 'Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'sign', theTerm.sign, ...
                    'formula', theTerm.expression, ...
                    'operators', theTerm.operators);
            elseif strcmp(theTerm.name, 'Implicit Source')
                theEquation = cfdAddTerm(theEquationName, theTerm.name, ...
                    'sign', theTerm.sign, ...
                    'formula', theTerm.expression);                
            end
        end
        
        % Update time
        theTerm = cfdGetTermInEquation(theEquationName, 'Transient');
        if isempty(theTerm)   
            cfdSetupMeshField(['rho_',theEquationName,'eq']);
            
            time = cfdGetTime;
            time.type = 'Steady';
            cfdSetTime(time);
        end
        
    end
end

% Setup Boundary Conditions
theMesh = cfdGetMesh;
theNumberOfBoundaries = theMesh.numberOfBoundaries;
for iBC = 1:theNumberOfBoundaries
    % Store boundary index
    theBC.index = iBC;
    
    % Get bc type
    boundaryField = theField.boundaryField{iBC};
    
    theBC.type = boundaryField.type;
    
    if strcmp(boundaryField.type, 'fixedValue') || strcmp(boundaryField.type, 'calculated') || strcmp(boundaryField.type, 'fixedGradient')        
        if strcmp(boundaryField.valueType, 'uniform')
            if strcmp(theFieldType, 'volVectorField')
                theBC.value = ['[',num2str(boundaryField.value(1)),';',num2str(boundaryField.value(2)),';',num2str(boundaryField.value(3)),']'];
            else
                theBC.value = num2str(boundaryField.value);
            end            
        else
            theBC.value = boundaryField.value;  
        end        
    elseif strcmp(boundaryField.type, 'zeroGradient') || strcmp(boundaryField.type, 'noSlip') || strcmp(boundaryField.type, 'empty') || strcmp(boundaryField.type, 'noSlip')
        if strcmp(theFieldType, 'volVectorField')
            theBC.value = '[0;0;0]';
        else
            theBC.value = '0';
        end
    end
    theEquation.bcs{iBC} = theBC;
end
    
% Store in data base
cfdSetModel(theEquation);





