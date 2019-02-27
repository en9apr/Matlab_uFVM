function list = readNonuniformList(key, fileDirectory, varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads the list from a nonuniform List
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
    
    if strcmp(key, 'internalField')
        C = textscan(tline, '%s nonuniform List%s', 1);
        
        if strcmp(C{1}{1}, key)
            % Read class (scalar or vector)
            fieldClass = C{2}{1};
            
            % Read list length
            C = textscan(fgetl(fileID), '%d', 1);
            listLength = C{1};
            
            fgetl(fileID);
            
            if strcmp(fieldClass, '<scalar>')
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '%f', 1);
                    list(i, 1) = C{1};
                end
            else
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '(%f %f %f)');
                    list(i, 1) = C{1};
                    list(i, 2) = C{2};
                    list(i, 3) = C{3};
                end
            end
        end
    else
        C = textscan(tline, '%s');
        if strcmp(C{1}{1}, key)
            
            while ~strcmp(C{1}{1}, 'value')
                C = textscan(fgetl(fileID), '%s nonuniform List%s', 1);
            end
            
            % Read class (scalar or vector)
            fieldClass = C{2}{1};            
            
            % Read list length
            C = textscan(fgetl(fileID), '%d', 1);
            listLength = C{1};
            
            fgetl(fileID);
            
            if strcmp(fieldClass, '<scalar>')
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '%f', 1);
                    list(i, 1) = C{1};
                end
            else
                for i=1:listLength
                    tline = fgetl(fileID);
                    C = textscan(tline, '(%f %f %f)');
                    list(i, 1) = C{1};
                    list(i, 2) = C{2};
                    list(i, 3) = C{3};
                end
            end
        end
    end
end

fclose(fileID);

