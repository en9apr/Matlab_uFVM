function theTerm = cfdSetupTerm(theTermName,theTermCoefficient,theTermVariable,theTermType,theTermSign)
%===================================================
%
% variableName ('')
% type (Residual) Deferred False
% sign (1) 
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
if(nargin ==1)
    theTermCoefficient = '';
    theTermVariable = ''; 
    theTermType = 'Residual';
    theTermSign = 1;
elseif(nargin == 2)
   theTermVariable = ''; 
    theTermType = 'Residual';
    theTermSign = 1;
elseif(nargin == 3)
    theTermType = 'Residual';
    theTermSign = 1;
elseif(nargin == 4)
    theTermSign = 1;
end


theTerm.name = theTermName;
theTerm.variableName = theTermVariable;
theTerm.coefficientName = theTermCoefficient;
theTerm.type = theTermType;
theTerm.sign = theTermSign;


if(strcmp(theTermName,'Source')==1)
   theTerm.Sb='0';
   theTerm.Sc='0';
   
elseif(strcmp(theTermName,'Convection')==1)
   theTerm.scheme = 'UPWIND';
    
elseif(strcmp(theTermName,'Transient')==1)
   theTerm.scheme = 'EULER';
    
end


end
