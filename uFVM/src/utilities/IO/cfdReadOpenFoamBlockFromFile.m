function blockData = cfdReadOpenFoamBlockFromFile(fileDirectory, blockName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads specific block from a FOAM file given the block
%   name and the file directory
%--------------------------------------------------------------------------
%
% Read File
fileID = fopen(fileDirectory, 'r');

% Skip Header
for i=1:16
    fgetl(fileID);
end

while(~feof(fileID))
    C = textscan(fileID, '%s', 1);    
    if(strcmp(C{1}, blockName))
        

    end    
end

error('\n%s\n',['"', blockName, '"',' block is not found in the file']);
















