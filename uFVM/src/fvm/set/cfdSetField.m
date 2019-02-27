function cfdSetField(theField)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theLocale = theField.locale;
theTimeStep = theField.timeStep;
theType = theField.type;

   theFieldIndex = 0;
    theNumberOfFields = length(Domain.(theTimeStep).(theLocale).fields);
    for iField=1:theNumberOfFields
       if(strcmp(theField.name,Domain.(theTimeStep).(theLocale).fields{iField}.name))
             Domain.(theTimeStep).(theLocale).fields{iField} = theField;
             theFieldIndex = iField;
       end
    end
    if(theFieldIndex==0) 
        Domain.(theTimeStep).(theLocale).fields{theNumberOfFields+1}=theField;
    end     
 

end