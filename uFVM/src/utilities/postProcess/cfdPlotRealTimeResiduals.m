function cfdPlotRealTimeResiduals(iter,iFigure)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function plots the real time residuals
%--------------------------------------------------------------------------

if(nargin==1)
    iFigure = 1;
end

theFigure = figure(iFigure);

axis tight
xlabel('Iterations');
if cfdIsTransient
    currentTime = cfdGetCurrentTime;
    title(['Time: ', num2str(currentTime), ' (s)']); 
end
%
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);
theLegend = {};
thePlot = [];
neqc = 0;

if iter==1
    theFigure.Name = 'Residuals';
    clf(iFigure);
    
    kk = 1;
    for iEquation=1:theNumberOfEquations
        theEquation = cfdGetModel(theEquationNames{iEquation});
        theEquationType = theEquation.type;
        theNumberOfComponents = 1;
        if strcmp(theEquationType,'Vector')
            theNumberOfComponents = 3;
        end
        for iComponent=1:theNumberOfComponents
            if(length(theEquation.terms) > 0)
                neqc = neqc + 1;
                if(neqc ==1) 
                    a='-xb';
                elseif(neqc ==2) 
                    a='-og';
                elseif(neqc ==3) 
                    a='-+r';
                elseif(neqc ==4) 
                    a='-<c';
                elseif(neqc ==5) 
                    a='-sm';
                elseif(neqc ==6) 
                    a='-dy';
                elseif(neqc ==7) 
                    a='-*k';
                end
                
                x = 1:iter;
                y = theEquation.residuals(x,iComponent)';
               
                h = semilogy(x, y, a, 'tag', ['h', num2str(kk)]);
                set(h, 'XDataSource', 'x');
                set(h, 'YDataSource', 'y');                
              
                thePlot = [thePlot, h];
                if strcmp(theEquationType,'Vector')
                    if iComponent==1
                        theEquationName = ['x-',theEquation.userName];
                    elseif iComponent==2
                        theEquationName = ['y-',theEquation.userName];
                    else
                        theEquationName = ['z-',theEquation.userName];
                    end                                        
                else
                    theEquationName = theEquation.userName;
                end                
                theLegend = {theLegend{:} theEquationName};
                hold on;
                
                kk = kk + 1;
            end
        end
    end
    legend(thePlot, theLegend);
    
    pause(1);
else
    kk = 1;
    for iEquation=1:theNumberOfEquations
        theEquation = cfdGetModel(theEquationNames{iEquation});
        theEquationType = theEquation.type;
        theNumberOfComponents = 1;
        if strcmp(theEquationType, 'Vector')
            theNumberOfComponents = 3;
        end
        for iComponent=1:theNumberOfComponents            
            x = 1:iter;
            y = theEquation.residuals(x,iComponent)';
            
            h = findobj('tag', ['h', num2str(kk)]);
                        
            refreshdata(h, 'caller');           
            
            kk = kk + 1;
        end        
    end
    drawnow;
    
end

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