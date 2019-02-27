function cfdSetupDomain
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes the database in which the variables will be
%   stored
%--------------------------------------------------------------------------
global Domain;
%-----------------------------------------
% Store Project Directory
%-----------------------------------------
Domain.projectDirectory = cd;
%-----------------------------------------
% Define Application Class
%-----------------------------------------
Domain.applicationClass = 'incompressible';
%-----------------------------------------
% Define MeshFields
%-----------------------------------------
Domain.Step0.Elements.fields = {};
Domain.Step0.Faces.fields = {};
Domain.Step0.Nodes.fields = {};
Domain.Step1.Elements.fields = {};
Domain.Step1.Faces.fields = {};
Domain.Step1.Nodes.fields = {};
%-----------------------------------------
% Define fluids, constants, properties 
%        equations, interfluidTerms
%-----------------------------------------
Domain.fluids = {};
Domain.interfluidTerms={};
Domain.constants = {};
Domain.fields = {};
%-----------------------------------------
% Define Fixed Element index
% <0 no element defined
%-----------------------------------------
Domain.iFixedElement = -1;
%-----------------------------------------
% Define Time Parameters
%-----------------------------------------
Domain.time = {};
%-----------------------------------------
% Define Basic Parameters
%-----------------------------------------
Domain.MAXResidual = 0.0001;
Domain.MAXIterations = 20;
Domain.currentIteration = 0;
%-----------------------------------------
% Define Basic Parameters
%-----------------------------------------
Domain.algorithm = 'SIMPLE';
%-----------------------------------------
% Define FOAM input
%-----------------------------------------
Domain.foam.controlDict = {};
Domain.foam.fvSchemes = {};
Domain.foam.fvSolution = {};
Domain.foam.transportProperties = {};
Domain.foam.thermophysicalProperties = {};
