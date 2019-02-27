function coefficients = cfdSetupCoefficients(theCConn,theCSize)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theNnumberOfElements = length(theCSize);

ac = zeros(theNnumberOfElements,1);
ac_old = zeros(theNnumberOfElements,1);
bc = zeros(theNnumberOfElements,1);

%dc & rc added to be used in ILU Solver
dc = zeros(theNnumberOfElements,1);
rc = zeros(theNnumberOfElements,1);

dphi = zeros(theNnumberOfElements,1);

anb = cell(theNnumberOfElements,1);

for iElement=1:theNnumberOfElements
    anb{iElement} = zeros(1,theCSize(iElement));
end

coefficients.ac = ac;
coefficients.ac_old = ac_old;
coefficients.bc = bc;
coefficients.anb = anb;

%dc & rc added to be used in ILU Solver
coefficients.dc = dc;
coefficients.rc = rc;

coefficients.dphi = dphi;

coefficients.cconn = theCConn;
coefficients.csize = theCSize;

coefficients.numberOfElements = theNnumberOfElements;
