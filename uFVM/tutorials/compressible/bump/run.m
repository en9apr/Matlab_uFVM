%--------------------------------------------------------------------------
%
%  written by the CFD Group @ AUB, 2017 
%  contact us at: cfd@aub.edu.lb
%==========================================================================
% Case Description:
%     In this test case a compressible flow against a bump at Mach 0.5 is
%     simulated
%--------------------------------------------------------------------------

% Setup Case
cfdSetupSolverClass('compressible');

% Read OpenFOAM Files
cfdReadOpenFoamFiles;

% Setup Time Settings
cfdSetupTime;

% Setup Equations
cfdDefineEquation('U', 'div(rho*U*U) =  laplacian(mu*U) - grad(p)'); % Momentum
cfdDefineEquation('p', 'div(rho*U) = 0'); % Continuity
cfdDefineEquation('T', 'div(rho*Cp*U*T) = laplacian(k*T)'); % Energy

% Initialize case
cfdInitializeCase;

% Run case
cfdRunCase;