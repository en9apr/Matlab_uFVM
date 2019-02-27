function value = cfdReadOpenFoamKeyValue(fileDirectory, key)
%===================================================

%  written by the CFD Group @ AUB, Fall 2017
%===================================================
%
% Read fvSchemes file
%
fpid = fopen(fileDirectory, 'r');

% Skip Header
for i=1:16
    tline = fgetl(fpid);
end

% Search for key name
splitted_str = strsplit(tline);
current_key = splitted_str{2};
while ~strcmp(current_key, key)
    tline = fgetl(fpid);
    splitted_str = strsplit(tline);
    if ~isempty(splitted_str{1})
        current_key = splitted_str{1};
    end
end

value = strjoin(splitted_str(2:end));
if strcmp(value(end),';')
    value(end) = [];
end










