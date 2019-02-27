function theInterfluidTerm = cfdSetupInterfluidTerm(theTermName,theTermModel,theFluid1Name,theFluid2Name)
%
%
theFluid1 = cfdGetFluidUsingName(theFluid1Name);
theFluid1Tag = theFluid1.tag;
theFluid1Index = theFluid1.index;

theFluid2 = cfdGetFluidUsingName(theFluid2Name);
theFluid2Tag = theFluid2.tag;
theFluid2Index = theFluid2.index;



theInterfluidTerm = fvmDefineInterfluidTerm(theTermName,theTermModel,theFluid1Index,theFluid2Index);

% ,theInterfluidTerm.iFluid1,theInterfluidTerm.iFluid2

% x component
theEquation1xName = ['Velx' theFluid1Tag];
theTerm1x=cfdAddTerm(theEquation1xName,'Drag');
theEquation2xName = ['Velx' theFluid2Tag];
theTerm2x=cfdAddTerm(theEquation2xName,'Drag');
% y component
theEquation1yName = ['Vely' theFluid1Tag];
theTerm1y=cfdAddTerm(theEquation1yName,'Drag');
theEquation2yName = ['Vely' theFluid2Tag];
theTerm2y=cfdAddTerm(theEquation2yName,'Drag');


theEquation1x = cfdGetModel(theEquation1xName);
theNumberOfTerms = length(theEquation1x.terms);

theEquation1x.terms{theNumberOfTerms}.index=theInterfluidTerm.index;
theEquation1x.terms{theNumberOfTerms}.iFluid1=theInterfluidTerm.iFluid1;
theEquation1x.terms{theNumberOfTerms}.iFluid2=theInterfluidTerm.iFluid2;

cfdSetModel(theEquation1x)



theEquation2x = cfdGetModel(theEquation2xName);
theNumberOfTerms = length(theEquation2x.terms);

theEquation2x.terms{theNumberOfTerms}.index=theInterfluidTerm.index;
theEquation2x.terms{theNumberOfTerms}.iFluid1=theInterfluidTerm.iFluid1;
theEquation2x.terms{theNumberOfTerms}.iFluid2=theInterfluidTerm.iFluid2;

cfdSetModel(theEquation2x)


theEquation1y = cfdGetModel(theEquation1yName);
theNumberOfTerms = length(theEquation1y.terms);

theEquation1y.terms{theNumberOfTerms}.index=theInterfluidTerm.index;
theEquation1y.terms{theNumberOfTerms}.iFluid1=theInterfluidTerm.iFluid1;
theEquation1y.terms{theNumberOfTerms}.iFluid2=theInterfluidTerm.iFluid2;

cfdSetModel(theEquation1y)

theEquation2y = cfdGetModel(theEquation2yName);
theNumberOfTerms = length(theEquation2y.terms);

theEquation2y.terms{theNumberOfTerms}.index=theInterfluidTerm.index;
theEquation2y.terms{theNumberOfTerms}.iFluid1=theInterfluidTerm.iFluid1;
theEquation2y.terms{theNumberOfTerms}.iFluid2=theInterfluidTerm.iFluid2;

cfdSetModel(theEquation2y)