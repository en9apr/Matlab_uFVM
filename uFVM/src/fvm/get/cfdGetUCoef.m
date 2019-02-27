function TM = cfdGetUCoef(iPatch)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates U coefficient
%--------------------------------------------------------------------------
global Domain;

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;

turbulenceProperties = Domain.foam.turbulenceProperties;
for iEntry=1:length(turbulenceProperties)
    if strcmp(turbulenceProperties{iEntry, 1}, 'turbulence')
        turbulence = Domain.foam.turbulenceProperties{iEntry, 2};
        break;
    end
end

walldist = [theMesh.faces(iBFaces).walldist]';
%
% Laminar Flow
%
if strcmp(turbulence,'off')
    theLaminarViscosityField = cfdGetMeshField('gamma_Ueq', 'Faces');
    visc_b = theLaminarViscosityField.phi(iBFaces);
    
    TM = visc_b./walldist;    
else
    % get turbulence model
    %
    % Turbulent Flow : k-epsilon model
    %
    for iEntry=1:length(turbulenceProperties)
        if strcmp(turbulenceProperties{iEntry, 1}, 'RASModel')
            RASModel = Domain.foam.turbulenceProperties{iEntry, 2};
            break;
        end
    end
    if strcmp(RASModel,'kEpsilon')        
        iOwners = [theMesh.faces(iBFaces).iOwner]';
        
        DensityField = cfdGetField('rho');
        density = DensityField.phi(iOwners);
        
        ViscosityField = cfdGetField('mu');
        visc = ViscosityField.phi(iOwners);
        
        TKEField = cfdGetField('TKE');
        tke = TKEField.phi(iOwners);
        
        Cmu = cfdGetConstant('Cmu');
        Kappa = cfdGetConstant('Kappa');
        Er = cfdGetConstant('E_roughness');
        
        
        u_star = Cmu^0.25 * tke.^0.5;
        y_plus = density.*u_star.*walldist./visc;
        u_plus = log10(y_plus*Er)/Kappa;
        
        TM = density.*u_star./u_plus;
        %
        % ERROR
        %
    else
        error('ERROR viscousModel not defined');
    end
end

end
