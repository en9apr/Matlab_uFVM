%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case the square cavity problem is considered with a
%     uniform velocity profile throughout the domain. The objective is to
%     investigate the convection schemes (the default now is set to first
%     order upwind).
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('basic');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('phi', 'div(rho*U*phi) = 0'); % Convection

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;