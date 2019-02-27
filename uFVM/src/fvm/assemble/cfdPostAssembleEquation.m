function cfdPostAssembleEquation(theEquationUserName,iComponent)
%===================================================
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
theEquationBaseName = cfdGetBaseName(theEquationUserName);

if(nargin==1)
    iComponent=1;
end
%
theCoefficients = cfdGetCoefficients;
%
theMesh = cfdGetMesh;
iElements = [1:theMesh.numberOfElements]';
%
ac  = theCoefficients.ac;
ac_old  = theCoefficients.ac_old;
anb = theCoefficients.anb;
bc  = theCoefficients.bc;

%---------------------------------------------------
% Momentum Equations
%---------------------------------------------------
if strcmp(theEquationBaseName,'U')
    
    theDUField = cfdGetMeshField(['DU' num2str(iComponent)]);
    theDUTField = cfdGetMeshField(['DUT' num2str(iComponent)]);
    
    %---------------------------------------------------------------------------
    %  Computations @ CELLS
    %---------------------------------------------------------------------------
    theAlgorithm = cfdGetAlgorithm;
    %
    % SIMPLE
    %
    if strcmp(theAlgorithm,'SIMPLE')
        volume = [theMesh.elements.volume]';
        
        theDUField.phi(iElements) = volume./ac;
        theDUTField.phi(iElements) = ac_old./ac;
    end
    %
    % SIMPLEC
    %
    if(strcmp(theAlgorithm,'SIMPLEC'))
        
        
    end
    %
    % SIMPLEST
    %
    if(strcmp(theAlgorithm,'SIMPLEST'))
        
        
    end
    %
    % SIMPLER
    %
    if(strcmp(theAlgorithm,'SIMPLER'))
        
        
    end
    %
    % PISO
    %
    if(strcmp(theAlgorithm,'PISO'))
        
        
    end
    
    %---------------------------------------------------------------------------
    %  Computations @ BOUNDARY
    %---------------------------------------------------------------------------
    iBElements = [theMesh.numberOfElements+1:theMesh.numberOfElements+theMesh.numberOfBFaces]';
    
    iBFaces= [theMesh.numberOfInteriorFaces+1:theMesh.numberOfFaces]';
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    
    theDUField.phi(iBElements) = theDUField.phi(iOwners);
    theDUTField.phi(iBElements) = theDUTField.phi(iOwners);
    
    %
    cfdSetMeshField(theDUField);
    cfdSetMeshField(theDUTField);      
    
end

theCoefficients.ac = ac;
theCoefficients.anb = anb;
theCoefficients.bc = bc;

cfdSetCoefficients(theCoefficients);

end
