%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case an air flow is simulated past a 2D cylinder
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('incompressible');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('U', 'div(U*U) = laplacian(nu*U) + div(nu*transp(grad(U))) - grad(p)'); % Momentum
cfdDefineEquation('p', 'div(U) = 0'); % Continuity

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;