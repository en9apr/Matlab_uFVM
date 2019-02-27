function theField = cfdGetField(theFullFieldName,theLocale,theTimeStep)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
global Domain;

if(nargin<3)
    theTimeStep = 'Step0';
end

if(nargin<2)
    theLocale = 'Elements';
end

theFieldName = cfdConvertName(theFullFieldName);

theField = {};
theNumberOfFields = length(Domain.(theTimeStep).(theLocale).fields);
for iField=1:theNumberOfFields
    theName = Domain.(theTimeStep).(theLocale).fields{iField}.name;
    if(strcmp(theFieldName,theName))
        theField = Domain.(theTimeStep).(theLocale).fields{iField};
        return;
    end
end
 