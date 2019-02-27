function cfdPrintResiduals(theEquationName,rmsResidual, maxResidual, lsResBefore, lsResAfter)

theNumberOfComponents = length(rmsResidual);

if theNumberOfComponents>1
    if length(theEquationName)>8
        theEquationName = theEquationName(1:8);
        RS = '';
        LS = '';
    elseif length(theEquationName)<8
        nSpaces = 8 - length(theEquationName);
        for nSpace=1:nSpaces
            spaces(nSpace)= ' ';
        end
        RS = spaces(1:floor(nSpaces/2));
        LS = spaces((floor(nSpaces/2)+1):end);
    else
        RS = '';
        LS = '';
    end
    
    comp = {'x','y','z'};
    
    for iComponent=1:theNumberOfComponents
        fprintf(['| ',RS,theEquationName,'-',comp{iComponent},LS, ...
            ' |  %.3E  |  %.3E  |    %.3E    |   %.3E   |\n'], ...
            rmsResidual(iComponent), maxResidual(iComponent), lsResBefore(iComponent), lsResAfter(iComponent));
    end
else
    
    if length(theEquationName)>10
        theEquationName = theEquationName(1:10);
        RS = '';
        LS = '';
    elseif length(theEquationName)<10
        nSpaces = 10 - length(theEquationName);
        for nSpace=1:nSpaces
            spaces(nSpace)= ' ';
        end
        RS = spaces(1:floor(nSpaces/2));
        LS = spaces((floor(nSpaces/2)+1):end);
    else
        RS = '';
        LS = '';
    end
    
    fprintf(['| ',RS,theEquationName,LS,' |  %.3E  |  %.3E  |    %.3E    |   %.3E   |\n'],rmsResidual,maxResidual,lsResBefore,lsResAfter);
end