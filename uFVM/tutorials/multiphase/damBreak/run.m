%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a dam break case is simulated
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('multiphase');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Define mixture properties
cfdSetupProperty('rho', 'recipe', 'Mixture');
cfdSetupProperty('mu', 'recipe', 'Mixture');

% Setup Equations
cfdDefineEquation('U', 'ddt(rho*U) + div(rho*U*U) = laplacian(mu*U) + div(mu*transp(grad(U))) - grad(p) + rho*g'); % Momentum
cfdDefineEquation('p', 'ddt(rho) + div(rho*U) = 0'); % Pressure
cfdDefineEquation('alpha1', 'ddt(rho1*alpha1) + div(rho1*U*alpha1) = 0'); % vof
cfdDefineEquation('alpha2', 'ddt(rho2*alpha2) + div(rho2*U*alpha2) = 0'); % vof

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;