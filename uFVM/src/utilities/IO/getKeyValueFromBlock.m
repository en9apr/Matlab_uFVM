function value = getKeyValueFromBlock(key, blockName, fileDirectory, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads they value of any entry (key) in a FOAM file given
%   the block name
%--------------------------------------------------------------------------

% Read File
fileID = fopen(fileDirectory, 'r');

if nargin==3
    skipHeader = true;
else
    skipHeader = varargin{1};
end

% Skip Header
if skipHeader
    for i=1:16
        fgetl(fileID);
    end
end

value = '';

% Check if required block is a subblock
nSubBlocks = length(strfind(blockName, '/'));
if nSubBlocks>0
    blockNames = textscan(blockName, '%s', nSubBlocks+1, 'delimiter', '/');
    blockNames = blockNames{1};
else
    blockNames = {blockName};
end

iBlock = 1;
while(~feof(fileID))
    % Read each line
    tline = fgetl(fileID);
    
    % Skip empty lines
    if isempty(tline)
        continue;
    end
    
    % Skip empty lines
    if isempty(strrep(tline, ' ', ''))
        continue;
    end
    
    % Skip commented lines
    if length(tline)>1
        if strcmp(tline(1:2), '//')
            continue;
        end
    end
    
    % Search for the block name
    C = textscan(tline, '%s', 1);
    if strcmp(C{1}{1}, blockNames{iBlock})
        if length(blockNames)~=iBlock
            iBlock = iBlock + 1;
            continue;
        end
    else
        continue;
    end
    
    % Skip to content
    tline = fgetl(fileID);
    tline = fgetl(fileID);
    
    % Collect value
    C = textscan(tline, '%s %[^\n]', 1);
    entry = C{1}{1};
    while ~strcmp(key, entry)
        tline = fgetl(fileID);
        
        if tline<0
            error('\n%s\n',[key,' block is not found in the file ', fileDirectory]);
        end
        
        % Skip empty lines
        if isempty(tline)
            continue;
        end
        
        % Skip empty lines
        if isempty(strrep(tline, ' ', ''))
            continue;
        end        
        
        % Skip commented lines
        if length(tline)>1
            if strcmp(tline(1:2), '//')
                continue;
            end
        end
        
        C = textscan(tline, '%s %[^\n]', 1);
        entry = C{1}{1};
    end
    semicolonLocation = strfind(C{2}{1}, ';');
    if ~isempty(semicolonLocation)
        value = C{2}{1}(1:semicolonLocation-1);
    else
        value = C{2}{1};
    end
    fclose(fileID);
    return;
end


