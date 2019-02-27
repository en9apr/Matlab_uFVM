function cfdInitializeCase
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes scalar fields
%--------------------------------------------------------------------------

% Print Action
fprintf('\nInitializing Fields ...\n\n\n');

global Domain;

% Check case. necessary
cfdCheckCase;

% Loop over equations and initialize each
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);
for iEquation = 1:theNumberOfEquations
    theEquationName = theEquationNames{iEquation};
    cfdInitializeEquation(theEquationName);
end

% Initialize Constant Fields. These fields are those that are defined in
% the 0  directory but are not included to be solved for (they dont have
% corresponding equations)
cfdInitializeConstantFields;

% Initialize Property Fields
thePropertyNames = cfdGetPropertyNames;
theNumberOfProperties = length(thePropertyNames);
for iProperty = 1:theNumberOfProperties
    thePropertyName = thePropertyNames{iProperty};
    cfdInitializeProperty(thePropertyName);
end

% Update Term Coefficient Fields (rho, psi and gamma)
for iEquation = 1:theNumberOfEquations
    theEquationName = theEquationNames{iEquation};
    
    % Exclude p equation from this step as it is not treated in a standard
    % way
    if strcmp(theEquationName, 'p')
        continue;
    end
    
    cfdUpdateTermCoefficientField(theEquationName);
    
    % If U equation, initialize mdot_f. If the current time step is not 0
    % directory, then mdot_f field is already in the foam data base as it
    % is collected from the mdot_f file in the time directory. It is then
    % created as a mesh field.
    if strcmp(theEquationName, 'U')
        if isfield(Domain.foam.fields, 'mdot_f')
            cfdCreateMeshField('mdot_f', Domain.foam.fields.mdot_f);
        else
            % If current time directory is 0, initialize the mdot_f field.
            % The foam field is to be created and stored in the foam data
            % base in order to use it to write the mdot_f in foam format to
            % the corresponding time directory
            cfdInitializeMdotField;
            foamField = cfdCreateFoamField('mdot_f', 'surfaceScalarField', '[1 0 -1 0 0 0 0]');
            Domain.foam.fields = setfield(Domain.foam.fields, 'mdot_f', foamField);
        end
    end
end

if cfdCheckIfFieldExists('U')
    theVelocityModel = cfdGetModel('U');
    if strcmp(theVelocityModel.class, 'Constant')
        cfdSetupMeshField('mdot_f','Faces');
        cfdCreateMdotField;
    end
end

%
% Start timing
%
Domain.startTime = double(tic);
Domain.currentTime = Domain.startTime;




