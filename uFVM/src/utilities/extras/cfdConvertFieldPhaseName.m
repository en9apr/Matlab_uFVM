function theFieldPhaseName = cfdConvertFieldPhaseName(theFieldPhaseName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function changes the field name from foam phase such that:
%   field.phase into fieldIndex
%--------------------------------------------------------------------------

global Domain;

if isfield(Domain.foam.transportProperties, 'phases')
    dotIndex = strfind(theFieldPhaseName, '.');
    if dotIndex
        phases = Domain.foam.transportProperties.phases;
        phaseNames = fieldnames(phases);
        for iPhase=1:length(phaseNames)
            if strcmp(theFieldPhaseName(dotIndex+1:end), phaseNames{iPhase})
                theFieldPhaseName(dotIndex:end) = '';
                theFieldPhaseName = [theFieldPhaseName, num2str(iPhase)];
                break;
            end
        end
    else
        phases = Domain.foam.transportProperties.phases;
        phaseNames = fieldnames(phases);        
        for iPhase=1:length(phaseNames)
            if iPhase<9
                if strcmp(theFieldPhaseName(end), num2str(iPhase))
                    theFieldPhaseName(end) = '';
                    theFieldPhaseName = [theFieldPhaseName, '.', phaseNames{iPhase}];
                    break;
                end
            end
        end        
    end
end