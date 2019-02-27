function thePropertyField = cfdGetThermoProperty(thePropertyName,theFluidTag,varargin)

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
numberOfBElements = theMesh.numberOfBElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theLocale = 'Elements';
iLocale = 1:numberOfElements+numberOfBElements;

if nargin>2
    theLocale = varargin{1};
    if strcmp(theLocale,'Faces')
        iLocale = 1:numberOfInteriorFaces;
    end
end
if nargin>3
    iLocale = varargin{2};
end

theDensityModel = cfdGetModel(['Density' theFluidTag]);
model = theDensityModel.model;

if strcmp(model,'Ideal Gas')
    if strcmp(thePropertyName,'drhodp')
        
        temperature = cfdGetMeshField(['Temperature' theFluidTag]);
        T = temperature.phi;
        gasConstant = cfdGetMeshField(['GasConstant' theFluidTag]);
        R = gasConstant.phi;
          
        drhodp = 1./(R.*T);
        if strcmp(theLocale,'Faces')
            thePropertyField = cfdInterpolateFromElementsToFaces('Average',thePropertyField);
        end
        thePropertyField = thePropertyField(iLocale);
        %
        %
        %
    elseif strcmp(thePropertyName,'Density')
        
        theFieldNames = cfdGetFieldNames;
        numberOfFields = cfdGetNumberOfFields;
        for iField=1:numberOfFields
            theFieldName = theFieldNames{iField};
            switch theFieldName
                case ['Temperature' theFluidTag]
                    temperature = cfdGetMeshField(theFieldName);
                    T = temperature.phi;
                case 'Pressure'
                    pressure = cfdGetMeshField(theFieldName);
                    P = pressure.phi;
                case ['GasConstant' theFluidTag]
                    gasConstant = cfdGetMeshField(theFieldName);
                    R = gasConstant.phi;
            end
        end
        
        thePropertyField = P./(R.*T);
        if strcmp(theLocale,'Faces')
            thePropertyField = cfdInterpolateFromElementsToFaces('Average',thePropertyField);
        end
        thePropertyField = thePropertyField(iLocale);
        %
        %
        %
    end
elseif strcmp(model,'Real Gas')
    
    if strcmp(thePropertyName,'drhodp')
        
        theFieldNames = cfdGetFieldNames;
        numberOfFields = cfdGetNumberOfFields;
        for iField=1:numberOfFields
            theFieldName = theFieldNames{iField};
            switch theFieldName
                case ['Temperature' theFluidTag]
                    temperature = cfdGetMeshField(theFieldName);
                    T = temperature.phi;
                case ['GasConstant' theFluidTag]
                    gasConstant = cfdGetMeshField(theFieldName);
                    R = gasConstant.phi;
            end
        end
        
        thePropertyField = 1./(R.*T);
        if strcmp(theLocale,'Faces')
            thePropertyField = cfdInterpolateFromElementsToFaces('Average',thePropertyField);
        end
        thePropertyField = thePropertyField(iLocale);
        %
        %
        %
    elseif strcmp(thePropertyName,'Density')
        
        theFieldNames = cfdGetFieldNames;
        numberOfFields = cfdGetNumberOfFields;
        for iField=1:numberOfFields
            theFieldName = theFieldNames{iField};
            switch theFieldName
                case ['Temperature' theFluidTag]
                    temperature = cfdGetMeshField(theFieldName);
                    T = temperature.phi;
                case 'Pressure'
                    pressure = cfdGetMeshField(theFieldName);
                    P = pressure.phi;
                case ['GasConstant' theFluidTag]
                    gasConstant = cfdGetMeshField(theFieldName);
                    R = gasConstant.phi;
            end
        end
        
        thePropertyField = P./(R.*T);
        if strcmp(theLocale,'Faces')
            thePropertyField = cfdInterpolateFromElementsToFaces('Average',thePropertyField);
        end
        thePropertyField = thePropertyField(iLocale);
        %
        %
        %
    end
    
end
