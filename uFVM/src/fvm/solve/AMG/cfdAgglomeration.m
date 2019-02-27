function [children, parents, CCONN] = cfdAgglomeration(anb,cconn,coarseningBy)

numberOfFineElements = size(anb,1);

minNumberOfAgglomeratedElements = coarseningBy;
maxNumberOfAgglomeratedElements = minNumberOfAgglomeratedElements + 3;
numberOfAgglomeratedElements = 1;
previousNumberOfAgglomeratedElements = 1;
%
%
ancestors = zeros(numberOfFineElements,1);
parents = zeros(numberOfFineElements,1);
%
for iElement=1:numberOfFineElements
    maxAnb(iElement) = max(-anb{iElement});
    parent(iElement) = 0;
end
%
beta = 0.5;
%
%
iParent = 1;
iSeed = 1;
iSeed_old = 1;
iSeed_previous = [];

while ~isempty(iSeed)    
    iGeneration = 1;  
    
    while numberOfAgglomeratedElements<minNumberOfAgglomeratedElements
        if iGeneration==1
            iElement = iSeed;
            ancestors(iSeed) = iSeed;
            children{iCoarseElement,1}(numberOfAgglomeratedElements) = iSeed;
            parents(iSeed) = iCoarseElement;
        end
        for nAncestor=1:length(iElement)
            iAncestor = iElement(nAncestor);
            for nNeighbour=1:length(cconn{iAncestor})
                iNeighbour = cconn{iAncestor}(nNeighbour);
                if parents(iNeighbour)==0
                    iS = find(cconn{iNeighbour}==iAncestor);
                    if anb{iAncestor}(nNeighbour)/anb_max>=beta && anb{iNeighbour}(iS)/anb_max>=beta
                        ancestors(iNeighbour) = iAncestor;
                        children{iCoarseElement,1}(numberOfAgglomeratedElements+1) = iNeighbour;
                        parents(iNeighbour) = iCoarseElement;
                        numberOfAgglomeratedElements = numberOfAgglomeratedElements + 1;
                    end
                end
            end
        end
        
        iGeneration = iGeneration + 1;
        
        if previousNumberOfAgglomeratedElements==numberOfAgglomeratedElements
            ancestors(children{iCoarseElement,1}) = 0;
            parents(children{iCoarseElement,1}) = 0;
            children(iCoarseElement,:) = [];
            iCoarseElement = iCoarseElement - 1;
            break;
        end
        previousNumberOfAgglomeratedElements = numberOfAgglomeratedElements;
        
        iElement = children{iCoarseElement,1};
        
    end
    
    if numberOfAgglomeratedElements>=maxNumberOfAgglomeratedElements
        iDroppedElements = maxNumberOfAgglomeratedElements+1:numberOfAgglomeratedElements;
        ancestors(children{iCoarseElement,1}(iDroppedElements)) = 0;
        parents(children{iCoarseElement,1}(iDroppedElements)) = 0;
        children{iCoarseElement,1}(iDroppedElements) = [];
    end
    
    numberOfAgglomeratedElements = 1;
    previousNumberOfAgglomeratedElements = 1;
    iCoarseElement = iCoarseElement + 1;
    
    iSeed_previous = [iSeed_previous iSeed];
    passed = false;
    for i=1:length(ancestors)
        if ancestors(i)==0
            iSeed = i;
            for ii=1:length(iSeed_previous)
                if iSeed==iSeed_previous(ii)
                    passed = true;
                    break;
                end
            end
            if passed
                continue;
            else
                break;
            end
        end
        passed = false;
    end
    
    if iSeed==iSeed_old
        break;
    end
    
    iSeed_old = iSeed;
end

for iElement=1:length(parents)
    if parents(iElement)==0
        iNeighbours = cconn{iElement};
        iParents = parents(iNeighbours);
        minSize = minNumberOfAgglomeratedElements;
        for i=1:length(iParents)
            if iParents(i)~=0
                currentMinSize = length(children{iParents(i)});
                if currentMinSize<=minSize
                    iMinSize = i;
                end
            end
        end
        children{iParents(iMinSize)}(end+1) = iElement;
        parents(iElement) = iParents(iMinSize);
    end
end

nCoarseNeighbour = 0;
available = false;

for nCoarseElement=1:size(children,1)
    for nChild=1:length(children{nCoarseElement,1})
        iChild = children{nCoarseElement}(nChild);
        iNeighbours = cconn{iChild};
        for nNeighbour=1:length(iNeighbours)
            iNeighbour = iNeighbours(nNeighbour);
            if parents(iNeighbour)~=nCoarseElement
                if nCoarseNeighbour>0
                    available = false;
                    for i=1:length(CCONN{nCoarseElement})
                        if CCONN{nCoarseElement}(i)==parents(iNeighbour)
                            available = true;
                            break;
                        end
                    end
                end
                if available==false
                    CCONN{nCoarseElement,1}(nCoarseNeighbour+1) = parents(iNeighbour);
                    nCoarseNeighbour = nCoarseNeighbour + 1;
                end
            end
        end
    end
    nCoarseNeighbour = 0;
    available = false;
end