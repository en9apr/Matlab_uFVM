function residual = cfdSORSolver
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
theCoefficients = cfdGetCoefficients;


ac  = theCoefficients.ac;
anb = theCoefficients.anb;
bc  = theCoefficients.bc;
numberOfElements = length(ac);

dphi = theCoefficients.dphi;

for iElement=1:numberOfElements
    cconn = theCoefficients.cconn{iElement};
    local_dphi = bc(iElement);
    for iLocalNeighbour = 1:length(cconn)
        iNeighbour = cconn(iLocalNeighbour);
        local_dphi = local_dphi - anb{iElement}(iLocalNeighbour)*dphi(iNeighbour);
    end
    dphi(iElement) = local_dphi/ac(iElement);
end

for iElement=numberOfElements:-1:1
    cconn = theCoefficients.cconn{iElement};
    local_dphi = bc(iElement);
    for iLocalNeighbour = 1:length(cconn)
        iNeighbour = cconn(iLocalNeighbour);
        local_dphi = local_dphi - anb{iElement}(iLocalNeighbour)*dphi(iNeighbour);
    end
    dphi(iElement) = local_dphi/ac(iElement);
end

theCoefficients.dphi = dphi;
cfdSetCoefficients(theCoefficients);

