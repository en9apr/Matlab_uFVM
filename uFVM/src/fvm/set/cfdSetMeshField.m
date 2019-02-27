function cfdSetMeshField(theMeshField,theTimeStep)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores the mesh field in the data base
%--------------------------------------------------------------------------

global Domain;

if(nargin<2)
    theTimeStep = 'Step0';
end

theLocale = theMeshField.locale;

theFieldIndex = 0;
theNumberOfFields = length(Domain.(theTimeStep).(theLocale).fields);
for iField=1:theNumberOfFields
    if(strcmp(theMeshField.name,Domain.(theTimeStep).(theLocale).fields{iField}.name))
        cfdApplyMeshFieldLimits(theMeshField);
        Domain.(theTimeStep).(theLocale).fields{iField} = theMeshField;
        theFieldIndex = iField;
    end
end
if(theFieldIndex==0)
    Domain.(theTimeStep).(theLocale).fields{theNumberOfFields+1} = theMeshField;
end


end