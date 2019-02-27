function cfdPlotCoarseElements(gridLevel,iCoarseElements)
%
%
%
theMesh = cfdGetMesh;
theNumberOfFineElements = theMesh.numberOfElements;
theCoarsestLevelCoefficients = cfdGetCoefficients(gridLevel);
theNumberOfCoarsestLevelElements = theCoarsestLevelCoefficients.numberOfElements;
childrenCell = cell(theNumberOfCoarsestLevelElements,1);

if(nargin==1)
    iCoarseElements = 1:theNumberOfCoarsestLevelElements;
end

for iLevel=1:gridLevel-1
    theCoefficients = cfdGetCoefficients(iLevel);
    theParents{iLevel,1} = theCoefficients.parents';
end

for iElement=1:theNumberOfFineElements
    iLevel = 1;
    iParent = iElement;
    while iLevel<gridLevel
        iParent = theParents{iLevel}(iParent);
        iLevel = iLevel + 1;
    end    
    childrenCell{iParent} = [childrenCell{iParent} iElement];    
end

%

axis off
axis equal

cfdPlotMesh('FaceAlpha', 0, 'EdgeAlpha', .1);
for iCoarseElement=iCoarseElements
    for iFineElement=childrenCell{iCoarseElement}        
        theElement = theMesh.elements(iFineElement);
        iFaces = theElement.iFaces;
        for iFace=iFaces
            theFace = theMesh.faces(iFace);
            theNodes = theMesh.nodes(theFace.iNodes);
            local_phi = iCoarseElement;
            XYZ = [theNodes.centroid]';
            h = patch(XYZ(:,1),XYZ(:,2),XYZ(:,3),local_phi,'FaceAlpha', 1, 'EdgeAlpha', .1);
        end
        
    end
end

view(30,40)

set(gca,'Clipping','off');
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