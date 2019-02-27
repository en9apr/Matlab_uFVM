function cfdClearResults
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function deletes the folders and files that have been generated
%   through the run
%--------------------------------------------------------------------------

timeSteps = cfdGetTimeSteps;
for timeStep=timeSteps'
   if timeStep~=0 
       if exist(num2str(timeStep), 'dir')
           [SUCCESS,MESSAGE,MESSAGEID] = rmdir(num2str(timeStep), 's');
           SUCCESS
       end
   end
end