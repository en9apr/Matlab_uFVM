function applicationClass = cfdGetApplicationClass
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function retrieves the class of the application in the database
%--------------------------------------------------------------------------
global Domain;

applicationClass = Domain.applicationClass;
