%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a blood flow in an artery is simulated
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('incompressible');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('U', 'div(rho*U*U) = laplacian(mu*U) + div(mu*transp(grad(U))) - grad(p) + rho*g'); % Momentum
cfdDefineEquation('p', 'div(U) = 0'); % Pressure

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;