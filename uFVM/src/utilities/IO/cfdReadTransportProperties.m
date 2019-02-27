function cfdReadTransportProperties
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads transport properties from
%   "constant/transportProperties" file
%--------------------------------------------------------------------------
global Domain;

% Check Class of Application
applicationClass = cfdGetApplicationClass;
if(strcmp(applicationClass, 'incompressible') || strcmp(applicationClass, 'basic'))
    
    fprintf('\nReading Transport Properties ...');
    
    projectDirectory = pwd;
    projectDirectory = strrep(projectDirectory, '\', '/');
    
    transportPropertiesFile = [projectDirectory, '/constant/transportProperties'];
    
    % Check if "transportProperties" exists
    if(~exist(transportPropertiesFile, 'file'))
        error('\n%s\n','"transportProperties" file is not found in "~foamDirectory/constant"');
    end
    
    % Collect constant transport properties from file
    [propertyNames, propertyDetails] = getKeyValue('', transportPropertiesFile);        
    
    % Store attributes in the foam data base and create property fields
    Domain.foam.transportProperties = {};
    for iProperty=1:length(propertyNames)
        duplicateIndex = strfind(propertyDetails{iProperty}, propertyNames{iProperty});
        if ~isempty(duplicateIndex)
            propertyDetails{iProperty}(duplicateIndex:duplicateIndex+length(propertyNames{iProperty})-1) = '';            
        end
        
        C = textscan(propertyDetails{iProperty}, '[%d %d %d %d %d %d %d] %f');
        
        % Store attirbutes
        Domain.foam.transportProperties = setfield(Domain.foam.transportProperties, ...
            propertyNames{iProperty}, {'dimensions', [C{1} C{2} C{3} C{4} C{5} C{6} C{7}]; ...
            'value', C{8}});
        
        % Create Corresponding mesh field
        propertyName = propertyNames{iProperty};
        property = getfield(Domain.foam.transportProperties, propertyName);
        propertyValue = num2str(property{2,2});
        cfdSetupProperty(propertyName, 'constant', propertyValue);
        
        fprintf(['\nReading ', propertyName,' ...']);
    end
    cfdSetupFluid('fluid');
elseif(strcmp(applicationClass, 'compressible'))
    
    fprintf('\nReading Themrophysical Properties ...');
    
    projectDirectory = pwd;
    projectDirectory = strrep(projectDirectory, '\', '/');
    
    thermophysicalPropertiesFile = [projectDirectory, '/constant/thermophysicalProperties'];
    
    % Check if "transportProperties" exists
    if(~exist(thermophysicalPropertiesFile, 'file'))
        error('\n%s\n','"thermophysicalProperties" file is not found in "~foamDirectory/constant"');
    end
    
    % Get thermo type
    transport = getKeyValueFromBlock('transport', 'thermoType', thermophysicalPropertiesFile);
    thermo = getKeyValueFromBlock('thermo', 'thermoType', thermophysicalPropertiesFile);
    equationOfState = getKeyValueFromBlock('equationOfState', 'thermoType', thermophysicalPropertiesFile);
    
    % Get molecular weight
    molWeight = getKeyValueFromBlock('molWeight', 'mixture/specie', thermophysicalPropertiesFile);
    Rbar = 8.314e3; % universal gas constant
    R = Rbar / str2double(molWeight);
    cfdSetupConstant('R', R);
    
    % Store transport model
    if strcmp(transport, 'const')
        mu = getKeyValueFromBlock('mu', 'mixture/transport', thermophysicalPropertiesFile);
        Domain.foam.thermophysicalProperties.transport{1, 1} = 'mu';
        Domain.foam.thermophysicalProperties.transport{1, 2} = mu;
        cfdSetupProperty('mu', 'constant', mu);
        
        Pr = getKeyValueFromBlock('Pr', 'mixture/transport', thermophysicalPropertiesFile);
        Domain.foam.thermophysicalProperties.transport{1, 2} = 'Pr';
        Domain.foam.thermophysicalProperties.transport{1, 2} = Pr;
        cfdSetupConstant('Pr', str2double(Pr));
    elseif strcmp(transport, 'sutherland')
        As = getKeyValueFromBlock('As', 'mixture/transport', thermophysicalPropertiesFile);
        Ts = getKeyValueFromBlock('Ts', 'mixture/transport', thermophysicalPropertiesFile);
        Domain.foam.thermophysicalProperties.transport{1, 1} = 'mu';
        Domain.foam.thermophysicalProperties.transport{1, 2} = 'As*sqrt(T)/(1 + Ts/T)';
        cfdSetupProperty('mu', 'model', [num2str(As), '*sqrt(T)/(1 + ', num2str(Ts), '/T)']);
    elseif strcmp(transport, 'polynomial')
        error('\n%s\n','polynomial transport model is not implemented');
    end
    
    % Store thermodynamics model
    if strcmp(thermo, 'hConst')
        Cp = getKeyValueFromBlock('Cp', 'mixture/transport', thermophysicalPropertiesFile);
        Domain.foam.thermophysicalProperties.thermo{1, 1} = 'Cp';
        Domain.foam.thermophysicalProperties.thermo{1, 2} = Cp;
        cfdSetupProperty('Cp', 'constant', Cp);
        
        Domain.foam.thermophysicalProperties.thermo{1, 1} = 'k';
        Domain.foam.thermophysicalProperties.thermo{1, 2} = 'Cp*mu/Pr';
        cfdSetupProperty('k', 'model', 'Cp.*mu./Pr');
    elseif strcmp(thermo, 'hPolynomial')
        error('\n%s\n','hPolynomial thermo model is not implemented');
    end
    
    % Store equationOfState model
    if strcmp(equationOfState, 'perfectGas')
        Domain.foam.thermophysicalProperties.equationOfState{1, 1} = 'rho';
        Domain.foam.thermophysicalProperties.equationOfState{1, 2} = 'p/(R*T)';
        cfdSetupProperty('rho', 'model', 'p./(R.*T)');
        
        % Create drhodp property, used in pressure equation assembly
        cfdSetupProperty('C_rho', 'model', '1./(R.*T)');
    elseif strcmp(transport, 'incompressiblePerfectGas')
        error('\n%s\n','incompressiblePerfectGas equationOfState model is not implemented');
    elseif strcmp(transport, 'Boussinesq')
        error('\n%s\n','Boussinesq equationOfState model is not implemented');
    end
    cfdSetupFluid('fluid', 'mode', 'compressible');
elseif(strcmp(applicationClass, 'multiphase'))
    
    fprintf('\nReading Transport Properties ...');
    
    projectDirectory = pwd;
    projectDirectory = strrep(projectDirectory, '\', '/');
    
    transportPropertiesFile = [projectDirectory, '/constant/transportProperties'];
    
    % Check if "transportProperties" exists
    if(~exist(transportPropertiesFile, 'file'))
        error('\n%s\n','"transportProperties" file is not found in "~foamDirectory/constant"');
    end
    
    Domain.foam.transportProperties = {};
    
    % Get phases
    [entry, value] = getKeyValue('phases', transportPropertiesFile);
    value = strtrim(value{1});
    value = value(2:end-1);
    phases = strsplit(value);
    
    for iPhase=1:length(phases)
        blockData = readBlockData(phases{iPhase}, transportPropertiesFile);
        for iBlockItem=1:length(blockData)
            item = blockData{iBlockItem, 1};
            if strcmp(item, 'transportModel')
                
                
            else
                propertyName = item;
                propertyDetails = blockData{iBlockItem, 2};
                C = textscan(propertyDetails, '[%d %d %d %d %d %d %d] %f');
                
                % Store attirbutes
                Domain.foam.transportProperties.phases.(phases{iPhase}).(propertyName).dimensions = [C{1} C{2} C{3} C{4} C{5} C{6} C{7}];
                Domain.foam.transportProperties.phases.(phases{iPhase}).(propertyName).value = C{8};
                
                % Create corresponding mesh field
                cfdSetupProperty([propertyName, num2str(iPhase)], 'constant', num2str(C{8}));
                
                fprintf(['\nReading ', propertyName,' of ', phases{iPhase},' ...']);
            end
        end
        cfdSetupFluid(phases{iPhase});
    end
    
elseif(strcmp(applicationClass, 'heatTransfer'))
    
    projectDirectory = pwd;
    projectDirectory = strrep(projectDirectory, '\', '/');
    
    transportPropertiesFile = [projectDirectory, '/constant/transportProperties'];
    thermophysicalPropertiesFile = [projectDirectory, '/constant/thermophysicalProperties'];
    
    % Check if "transportProperties" exists
    if(~exist(transportPropertiesFile, 'file') && ~exist(thermophysicalPropertiesFile, 'file'))
        error('\n%s\n','"transportProperties" and/or "thermophysicalPropertiesFile" file are not found in "~foamDirectory/constant"');
    end
    
    if (exist(transportPropertiesFile, 'file')) 
        
        fprintf('\nReading Transport Properties ...');
        
        % Collect constant transport properties from file
        [propertyNames, propertyDetails] = getKeyValue('', transportPropertiesFile);
        
        % Store attributes in the foam data base and create property fields
        Domain.foam.transportProperties = {};
        for iProperty=1:length(propertyNames)
            if strcmp(propertyNames{iProperty}, 'transportModel')
                continue;
            end
            
            C = textscan(propertyDetails{iProperty}, '[%d %d %d %d %d %d %d] %f');
            
            % Store attirbutes
            Domain.foam.transportProperties = setfield(Domain.foam.transportProperties, ...
                propertyNames{iProperty}, {'dimensions', [C{1} C{2} C{3} C{4} C{5} C{6} C{7}]; ...
                'value', C{8}});
            
            % Create Corresponding mesh field
            propertyName = propertyNames{iProperty};
            property = getfield(Domain.foam.transportProperties, propertyName);
            propertyValue = num2str(property{2,2});
            
            if strcmp(propertyName, 'TRef') || ...
                    strcmp(propertyName, 'beta') || ...
                    strcmp(propertyName, 'Pr')
                constant.name = propertyName;
                constant.value = str2double(propertyValue);
                cfdSetConstant(constant);
            else
                cfdSetupProperty(propertyName, 'constant', propertyValue);
            end
            fprintf(['\nReading ', propertyName,' ...']);
        end
        cfdSetupFluid('fluid');
    end
    
    % Read thermophysical properties file if exists
    if exist(thermophysicalPropertiesFile, 'file')
        
        fprintf('\nReading Themrophysical Properties ...');
        
        % Get thermo type
        transport = getKeyValueFromBlock('transport', 'thermoType', thermophysicalPropertiesFile);
        thermo = getKeyValueFromBlock('thermo', 'thermoType', thermophysicalPropertiesFile);
        equationOfState = getKeyValueFromBlock('equationOfState', 'thermoType', thermophysicalPropertiesFile);
        
        % Get molecular weight
        molWeight = getKeyValueFromBlock('molWeight', 'mixture/specie', thermophysicalPropertiesFile);
        Rbar = 8.314e3; % universal gas constant
        R = Rbar / str2double(molWeight);
        cfdSetupConstant('R', R);
        
        % Store transport model
        if strcmp(transport, 'const')
            mu = getKeyValueFromBlock('mu', 'mixture/transport', thermophysicalPropertiesFile);
            Domain.foam.thermophysicalProperties.transport{1, 1} = 'mu';
            Domain.foam.thermophysicalProperties.transport{1, 2} = mu;
            cfdSetupProperty('mu', 'constant', mu);
            
            Pr = getKeyValueFromBlock('Pr', 'mixture/transport', thermophysicalPropertiesFile);
            Domain.foam.thermophysicalProperties.transport{1, 2} = 'Pr';
            Domain.foam.thermophysicalProperties.transport{1, 2} = Pr;
            cfdSetupConstant('Pr', str2double(Pr));
        elseif strcmp(transport, 'sutherland')
            As = getKeyValueFromBlock('As', 'mixture/transport', thermophysicalPropertiesFile);
            Ts = getKeyValueFromBlock('Ts', 'mixture/transport', thermophysicalPropertiesFile);
            Domain.foam.thermophysicalProperties.transport{1, 1} = 'mu';
            Domain.foam.thermophysicalProperties.transport{1, 2} = 'As*sqrt(T)/(1 + Ts/T)';
            cfdSetupProperty('mu', 'model', [num2str(As), '*sqrt(T)/(1 + ', num2str(Ts), '/T)']);
        elseif strcmp(transport, 'polynomial')
            error('\n%s\n','polynomial transport model is not implemented');
        end
        
        % Store thermodynamics model
        if strcmp(thermo, 'hConst')
            Cp = getKeyValueFromBlock('Cp', 'mixture/transport', thermophysicalPropertiesFile);
            Domain.foam.thermophysicalProperties.thermo{1, 1} = 'Cp';
            Domain.foam.thermophysicalProperties.thermo{1, 2} = Cp;
            cfdSetupProperty('Cp', 'constant', Cp);
            
            Domain.foam.thermophysicalProperties.thermo{1, 1} = 'k';
            Domain.foam.thermophysicalProperties.thermo{1, 2} = 'Cp*mu/Pr';
            cfdSetupProperty('k', 'model', 'Cp.*mu./Pr');
        elseif strcmp(thermo, 'hPolynomial')
            error('\n%s\n','hPolynomial thermo model is not implemented');
        end
        
        % Store equationOfState model
        if strcmp(equationOfState, 'perfectGas')
            Domain.foam.thermophysicalProperties.equationOfState{1, 1} = 'rho';
            Domain.foam.thermophysicalProperties.equationOfState{1, 2} = 'p/(R*T)';
            cfdSetupProperty('rho', 'model', 'p./(R.*T)');
            
            % Create drhodp property, used in pressure equation assembly
            cfdSetupProperty('C_rho', 'model', '1./(R.*T)');
        elseif strcmp(transport, 'incompressiblePerfectGas')
            error('\n%s\n','incompressiblePerfectGas equationOfState model is not implemented');
        elseif strcmp(transport, 'Boussinesq')
            error('\n%s\n','Boussinesq equationOfState model is not implemented');
        end
        cfdSetupFluid('fluid', 'mode', 'compressible');                
    end    
end
