%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a hot room is simulated with boussinesq
%     approximation
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('heatTransfer');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Define new properties
cfdSetupProperty('alpha', 'model', 'nu/Pr');

% Setup Equations
cfdDefineEquation('U', 'div(U*U) = laplacian(nu*U) + div(nu*transp(grad(U))) - grad(p) - g*beta*(T - TRef)'); % Momentum
cfdDefineEquation('p', 'div(U) = 0'); % Continuity
cfdDefineEquation('T', 'div(U*T) = laplacian(alpha*T)'); % Energy

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;