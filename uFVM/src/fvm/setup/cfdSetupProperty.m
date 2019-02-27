function cfdSetupProperty(theUserName,varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function set up properties
%--------------------------------------------------------------------------
%
global Domain;

if(isfield(Domain,'fields')==0)
    Domain.fields = {};
end

theNumberOfSubTerms = length(varargin);

thePropertyName = cfdConvertName(theUserName);

theProperty.userName = theUserName;
theProperty.name = thePropertyName;
theProperty.class = 'Property';
theProperty.type = 'Scalar';
theProperty.urf = 1.0;
theProperty.calcGradient = false;
theProperty.gradientType = '';


for iSubTerm = 1:2:theNumberOfSubTerms
   theProperty = setfield(theProperty,varargin{iSubTerm},varargin{iSubTerm+1});
end

if(theProperty.calcGradient)
    theProperty.phiGradient = [];
end

%
% Setup Associate MeshField
%
cfdSetupMeshField(theUserName,'Elements',theProperty.type);

cfdSetModel(theProperty);