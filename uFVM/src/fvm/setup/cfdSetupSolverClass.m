function cfdSetupSolverClass(applicationClass)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% cfdSetupCase   Setup new uFVM case.
%     cfdSetupCase(applicationClass) defines a new uFVM case of class by
%     creating the data base.
%
%     Examples:
%
%     % Create case of class incompressible:
%     cfdSetupSolverClass('incompressible')
%--------------------------------------------------------------------------

% Check if application class falls within the defined ones
if strcmp(applicationClass, 'incompressible') || ...
        strcmp(applicationClass, 'compressible') || ...
        strcmp(applicationClass, 'basic') || ...
        strcmp(applicationClass, 'multiphase') || ...
        strcmp(applicationClass, 'heatTransfer')
    fprintf(['\nRunning ', applicationClass,' ...\n']);
else
    error('\n%s\n',['"',applicationClass,'"',' application is not defined.', ...
        ' The following applications are applicable:'], 'basic', ...
        'incompressible', ...
        'compressible', ...
        'multiphase', ...
        'heatTransfer');
end

% Create Domain, a global variable, which is the data base of all case
cfdSetupDomain;

% Setup application class
cfdSetApplicationClass(applicationClass);