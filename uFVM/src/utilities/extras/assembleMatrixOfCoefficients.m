function [A, b] = assembleMatrixOfCoefficients

theCoefficients = cfdGetCoefficients;
ac = theCoefficients.ac;
bc = theCoefficients.bc;
anb = theCoefficients.anb;
cconn = theCoefficients.cconn;
numberOfElements = theCoefficients.numberOfElements;
A = zeros(numberOfElements);

for iElement=1:numberOfElements
   A(iElement, iElement) = ac(iElement);
   A(iElement, cconn{iElement}) = anb{iElement};
end
b = bc;

