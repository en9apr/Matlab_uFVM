function residual = cfdComputeResidualsArray(theCoefficients)

ac = theCoefficients.ac;
anb = theCoefficients.anb;
bc = theCoefficients.bc;
cconn = theCoefficients.cconn;
phi = theCoefficients.dphi;

numberOfElements = length(ac);
residual = zeros(size(ac));
for iElement=1:numberOfElements    
    residual(iElement) = bc(iElement) - ac(iElement)*phi(iElement);    
    for nNeighbour=1:length(cconn{iElement})
        iNeighbour = cconn{iElement}(nNeighbour);
        residual(iElement) = residual(iElement) - anb{iElement}(nNeighbour)*phi(iNeighbour);
    end    
end