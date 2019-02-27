function cfdSetApplicationClass(varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores the class of the application in the database
%--------------------------------------------------------------------------
global Domain;
if nargin==0    
    list = dir('..');
    [directory, foldername] = fileparts(list(1).folder);
    Domain.applicationClass = foldername;
else
    Domain.applicationClass = varargin{1};
end
