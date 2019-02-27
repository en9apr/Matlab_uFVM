%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a water flow in an elbow is simulated
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('incompressible');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('U', 'div(rho*U*U) = laplacian(mu*U) + div(mu*transp(grad(U))) - grad(p)'); % Momentum
cfdDefineEquation('p', 'div(U) = 0'); % Continuity

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;