function cfdPlotResiduals(iFigure)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

if(nargin==0)
    iFigure = 1;
end
theFigure = figure(iFigure);
theFigure.Name = 'Residuals';
clf(iFigure);

%
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);
theLegend = {};
thePlot = [];
neqc=0;
for iEquation=1:theNumberOfEquations
    theEquation = cfdGetModel(theEquationNames{iEquation});
    theEquationType = theEquation.type;
    theNumberOfComponents = 1;
    if(strcmp(theEquationType,'Vector'))
        theNumberOfComponents=3;
    end;
    for iComponent=1:theNumberOfComponents
        if(length(theEquation.terms) > 0)
            neqc = neqc+1;
            if(neqc ==1) a='-xb';end
            if(neqc ==2) a='-og';end
            if(neqc ==3) a='-+r';end
            if(neqc ==4) a='-<c';end
            if(neqc ==5) a='-sm';end
            if(neqc ==6) a='-dy';end
            if(neqc ==7) a='-*k';end
            theLength=length(theEquation.residuals(:,iComponent));
            h1=semilogy([1:theLength],theEquation.residuals(:,iComponent)',a);
            thePlot=[thePlot,h1];
            theEquationName = [theEquation.userName '(' num2str(iComponent) ')'];
            theLegend = {theLegend{:} theEquationName };
            hold on;
        end
    end;
end
legend(thePlot,theLegend);

axis tight
%ylim([1e-6 1e5])

end

function rgb = get_rgb(colour)
switch lower(colour),
    case {'y', 'yellow' }, rgb = [1, 1, 0];
    case {'m', 'magenta'}, rgb = [1, 0, 1];
    case {'c', 'cyan'   }, rgb = [0, 1, 1];
    case {'r', 'red'    }, rgb = [1, 0, 0];
    case {'g', 'green'  }, rgb = [0, 1, 0];
    case {'b', 'blue'   }, rgb = [0, 0, 1];
    case {'w', 'white'  }, rgb = [1, 1, 1];
    case {'k', 'black'  }, rgb = [0, 0, 0];
    otherwise            , rgb = [0, 0, 1]; % Unknown colour -> 'blue'.
end
end