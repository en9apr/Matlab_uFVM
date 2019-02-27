function cfdCreateConstantField(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function creates constant fields
%--------------------------------------------------------------------------

global Domain;

fieldNames = fieldnames(Domain.foam.fields);

for iField=1:length(fieldNames)
    if strcmp(theFieldName, fieldNames{iField})
        % get field
        theField = getfield(Domain.foam.fields, theFieldName);
        
        % Get formula
        theFieldIC = theField.internalField;
        if strcmp(theFieldIC.valueType, 'uniform')
            if strcmp(theField.class, 'volVectorField')
                ic = value;
            else
                C = textscan(theFieldIC, '%s %f');
                theEquation.ic = num2str(C{2});
            end
        else
            ic = theFieldIC.value;
        end
        
        % Evaluate at mesh
        if strcmp(theField.class, 'volScalarField')
            cfdSetupMeshField(theFieldName,'Elements','Scalar');
        else
            cfdSetupMeshField(theFieldName,'Elements','Vector');
        end
        phi = ic;
        
        theMeshField = cfdGetMeshField(theFieldName);
        theMeshField.phi = phi;
        cfdSetMeshField(theMeshField);
        
        % Evaluate at boundary
        theMesh = cfdGetMesh;
        theNumberOfBoundaries = theMesh.numberOfBoundaries;
        for iBC = 1:theNumberOfBoundaries
            % Store boundary index
            theBC.index = iBC;
            
            % Get bc type
            theBoundaryField = theField.boundaryField{iBC};
            
            % Initialize the boundary with a fixedValue type and uniform value
            theBC.type = 'fixedValue';
            if strcmp(theField.class, 'Scalar')
                theBC.value = '0';
            elseif strcmp(theField.class, 'Vector')
                theBC.value = '[0;0;0]';
            end
            
            % Retrieve from data base
            for iEntry=1:size(theBoundaryField, 1)
                if strcmp(theBoundaryField{iEntry, 1}, 'type')
                    theBC.type = theBoundaryField{iEntry, 2};
                elseif strcmp(theBoundaryField{iEntry, 1}, 'value')
                    if strfind(theBoundaryField{iEntry, 2}, 'uniform')
                        if strcmp(theField.class, 'volVectorField')
                            op = strfind(theBoundaryField{iEntry, 2}, '(');
                            cp = strfind(theBoundaryField{iEntry, 2}, ')');
                            C = textscan(theBoundaryField{iEntry, 2}(op+1:cp-1), '%f %f %f');
                            theBC.value = ['[',num2str(C{1}),';',num2str(C{2}),';',num2str(C{3}),']'];
                        else
                            C = textscan(theBoundaryField{iEntry, 2}, '%s %f');
                            theBC.value = num2str(C{2});
                        end
                    else
                        theBC.value = theBoundaryField{iEntry, 2};
                    end
                else
                    theBC = setfield(theBC, theBoundaryField{iEntry, 1}, theBoundaryField{iEntry, 2});
                end
            end
            
            bcs{iBC} = theBC;
        end
        
        if strcmp(theField.class, 'Scalar')
            cfdEvaluateConstantFieldBoundaryConditions(theFieldName, bcs, 'Scalar');
        elseif strcmp(theField.class, 'Vector')
            cfdEvaluateConstantFieldBoundaryConditions(theFieldName, bcs, 'Vector');
        end
                
        break;
    end
end
