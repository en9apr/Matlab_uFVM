function [AC, ANB] = cfdGetCoarseLevelCoefficients(children, parents, CCONN, ac, anb, cconn)

numberOfCoarseElements = size(children,1);
AC = zeros(numberOfCoarseElements,1);

for nCoarseElement=1:numberOfCoarseElements
    ANB{nCoarseElement,1} = CCONN{nCoarseElement}*0;
    for nFineElement=1:length(children{nCoarseElement})
        iFineElement = children{nCoarseElement}(nFineElement);
        
        AC(nCoarseElement) = AC(nCoarseElement) + ac(iFineElement);
        
        for nNeighbourFineElement=1:length(cconn{iFineElement})
            iNeighbourFineElement = cconn{iFineElement}(nNeighbourFineElement);
            if parents(iNeighbourFineElement)==nCoarseElement
                nF = find(cconn{iNeighbourFineElement}==iFineElement);
                
                AC(nCoarseElement) = AC(nCoarseElement) + ...
                    anb{iFineElement}(nNeighbourFineElement) + ...
                    anb{iNeighbourFineElement}(nF);
            else
                nFF = find(CCONN{nCoarseElement}==parents(iNeighbourFineElement));
                ANB{nCoarseElement}(nFF) = ANB{nCoarseElement}(nFF) + anb{iFineElement}(nNeighbourFineElement);
            end            
        end
    end        
end