function cfdSetupTime
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets up the time settings
%--------------------------------------------------------------------------

global Domain;
controlDict = Domain.foam.controlDict;
ddtSchemes = Domain.foam.fvSchemes.ddtSchemes;

for iEntry=1:size(controlDict,1)
    switch controlDict{iEntry, 1}
        case 'startFrom'
            time.startFrom = controlDict{iEntry, 2};
        case 'startTime'
            if strcmp(time.startFrom, 'firstTime')
                time.startTime = 0;
            elseif strcmp(time.startFrom, 'latestTime')
                timeDirectories = cfdGetTimeSteps;
                time.startTime = max(timeDirectories);
            elseif strcmp(time.startFrom, 'startTime')
                time.startTime = controlDict{iEntry, 2};
            end                                       
        case 'stopAt'
            time.stopAt = controlDict{iEntry, 2};            
        case 'endTime'
            time.endTime = controlDict{iEntry, 2};
        case 'deltaT'
            time.deltaT = controlDict{iEntry, 2};           
    end
end

time.fdt = 1e9; % false transience. Default

% Check ddtSchemes dictionary
if strcmp(ddtSchemes{2}, 'steadyState')
    time.type = 'Steady';
else
    time.type = 'Transient';
end

% Store time scheme in time data base
cfdSetTime(time);