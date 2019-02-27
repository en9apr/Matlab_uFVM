function [RMSResidual, MAXResidual] = cfdComputeNormalizedResidual(theEquationUserName, iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function computes normalized residuals
%--------------------------------------------------------------------------

theEquationBaseName = cfdGetBaseName(theEquationUserName);
theEquationField = cfdGetMeshField(theEquationUserName);

phiRho = [];

theCoefficients = cfdGetCoefficients;
ac = theCoefficients.ac;
anb = theCoefficients.anb;
bc = theCoefficients.bc;
cconn = theCoefficients.cconn;

theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
theElements = theMesh.elements;

% Residual for other Fields
if strcmp(theEquationBaseName,'p')
    effDiv = cfdComputeEffectiveDivergence;
    MAXResidual = max(effDiv);
    RMSResidual = sqrt(sum(effDiv.^2)/length(effDiv));
else
    
    %     % Loop over elements and calculate residual at each element
    %     Rc = abs(bc);
    %
    %     % Residuals. Calculate for convenience. Otherwise, they are not used
    %     Rc_max = max(Rc);
    %     Rc_rms = sqrt(sum(Rc.^2)/theNumberOfElements);
    %
    %     % Get phi scale from data base.
    %     % phi_scale = max(abs(phi)). And if pgi is zero, phi_scale is set to 1
    phiScale = cfdGetScale(theEquationUserName);
    %
    %     % Normalized Residuals
    %     Rc_scaled = Rc / (max(abs(ac))*phiScale);
    %     Rc_max_scaled = max(Rc_scaled);
    %     Rc_rms_scaled = sqrt(sum(Rc_scaled.^2)/theNumberOfElements);
    %
    %     MAXResidual = Rc_max_scaled;
    %     RMSResidual = Rc_rms_scaled;
    
    
    % Another approach which takes the transient term into consideration
    %
    theTransientTerm = cfdGetTermInEquation(theEquationUserName,'Transient');
    isTransient = cfdIsTransient && (~isempty(theTransientTerm));
    %
    if isTransient
        theRhoField = cfdGetMeshField(['rho_',theEquationBaseName,'eq']);
        phiRho = theRhoField.phi;
        dt = cfdGetDt;
    end
    
    theResidualSquared = 0.;
    theMaxResidual = 0;
    for iElement=1:theNumberOfElements
        volume = theElements(iElement).volume;
        local_ac = ac(iElement);
        if isTransient
            at = volume*phiRho(iElement)/dt;
            local_ac = local_ac - at;
            if(local_ac < 1e-6*at)
                local_ac = at;
            end
        end
        local_residual = bc(iElement);
        if ac(iElement) == 0
            disp(iElement);
            disp(theEquationField.name);
            error('****** Divide by Zero *******');
        end
        local_residual = local_residual/(local_ac*phiScale);
        theMaxResidual = max(theMaxResidual,abs(local_residual));
        theResidualSquared = theResidualSquared + local_residual*local_residual;
    end
    %
    MAXResidual = theMaxResidual;
    RMSResidual = sqrt(theResidualSquared/theNumberOfElements);
end

end