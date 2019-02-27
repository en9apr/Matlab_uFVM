function theNumberOfAMGLevels = cfdGetNumberOfAMGLevels
%
% Get number of AMG Levels for the size of Coefficients
%
global Domain;

theNumberOfAMGLevels = length(Domain.coefficients);