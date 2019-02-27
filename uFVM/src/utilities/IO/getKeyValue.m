function [entry, value] = getKeyValue(key, fileDirectory, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads they value of any entry (key) in a FOAM file
%--------------------------------------------------------------------------
%
% Read File
fileID = fopen(fileDirectory, 'r');

if nargin==2
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

if ~iscell(key)
    key = {key};
end

value = '';

matched_keys = {};
iMatched_keys = 1;
iEntry = 1;

while(~feof(fileID))
    % Read each line
    tline = fgetl(fileID);
    
    % Skip empty lines
    if isempty(tline)
        continue;
    end  
        
    % Skip commented lines
    if length(tline)>1
        if strcmp(tline(1:2), '//')
            continue;
        end
    end
    
    if strcmp(key{1}, '')
        C = textscan(tline, '%s %[^\n]');
        entry{iEntry} = C{1}{1};
        value{iEntry} = C{2}{1}(1:end-1);
        iEntry = iEntry + 1;
    else
        % Read first token
        C = textscan(tline, '%s', 1);
        
        for iKey=1:length(key)
            % If key matches first token, get the value
            if strcmp(C{1}{1}, key{iKey})
                matched_keys{iMatched_keys} = key{iKey};
                iMatched_keys = iMatched_keys + 1;
                
                C = textscan(tline, '%*s %[^\n]');
                
                % return empty string if no value exists
                if isempty(C{1})
                    value{iKey} = '';
                    entry{iKey} = key{iKey};
                    continue;
                end
                
                value{iKey} = C{1}{1};
                entry{iKey} = key{iKey};
                
                % Remove semi-colon if exists at the end of the value
                if strcmp(value{iKey}(end), ';')
                    value{iKey}(end) = [];
                end
            end
        end
    end
end

if ~strcmp(key{1}, '')
    if length(matched_keys)~=length(key)
        iUnmatched_keys = 1;
        for iKey=1:length(key)
            if ~any(ismember(matched_keys, key{iKey}))
                unmatched_keys{iUnmatched_keys} = key{iKey};
                iUnmatched_keys = iUnmatched_keys + 1;
            end
        end
        
        if length(unmatched_keys)>1
            error('\n%s\n',[strjoin(unmatched_keys(1:end-1), ', '), ' and ', unmatched_keys{end},' entries are not found in the file "', fileDirectory,'"']);
        else
            error('\n%s\n',[unmatched_keys{1},' entry is not found in the file "', fileDirectory,'"']);
        end
    end
end

fclose(fileID);
