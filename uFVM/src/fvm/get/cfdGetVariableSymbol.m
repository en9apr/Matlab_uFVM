function theSymbol = cfdGetVariableSymbol(theName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2017
%===================================================

cPos = strfind(theName,':');
if(isempty(cPos))
    theBaseName = theName;
elseif(length(cPos)==1)    
    theBaseName = theName(1:cPos-1);
elseif(length(cPos)==2)    
    theBaseName = theName(1:cPos(1)-1);
end

if(strcmp(theBaseName, 'Density'))
    theSymbol = 'rho';
elseif (strcmp(theBaseName, 'Viscosity'))
    theSymbol = 'mu';
elseif (strcmp(theBaseName, 'SpecificHeat'))
    theSymbol = 'Cp';
elseif (strcmp(theBaseName, 'SpecificHeatRatio'))     
    theSymbol = 'gamma';
elseif (strcmp(theBaseName, 'ThermalConductivity'))
    theSymbol = 'k';
elseif (strcmp(theBaseName, 'Velocity'))    
    theSymbol = 'U';
elseif (strcmp(theBaseName, 'Pressure'))    
    theSymbol = 'p';
elseif (strcmp(theBaseName, 'Temperature'))
    theSymbol = 'T';
elseif (strcmp(theBaseName, 'GasConstant'))        
    theSymbol = 'R';
else
    theSymbol = theBaseName;
end

% If more than one fluid is available, add an index at the end of every
% tag name
theNumberOfFluids = cfdGetNumberOfFluids;
theFluidName = cfdGetFluidName(theName);
if theNumberOfFluids>1 && ~isempty(theFluidName)
    theSymbol = [theSymbol, '_', theFluidName(1)];
end
