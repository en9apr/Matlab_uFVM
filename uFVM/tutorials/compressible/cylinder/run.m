%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case an air flow is simulated past a 2D cylinder
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('compressible');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('U', 'div(rho*U*U) = laplacian(mu*U) + div(mu*transp(grad(U))) - grad(p) - 2/3*grad(mu*div(U))'); % Momentum
cfdDefineEquation('p', 'div(rho*U) = 0'); % Continuity
cfdDefineEquation('T', 'div(rho*Cp*U*T) = laplacian(k*T) + DDt(p)'); % Energy

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;