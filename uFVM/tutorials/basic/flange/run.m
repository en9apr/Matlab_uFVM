%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a temperature diffusion of water flow in a square
%     cavity is simulated
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('basic');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('T', 'laplacian(DT*T) = 0'); % Heat Diffusion

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;