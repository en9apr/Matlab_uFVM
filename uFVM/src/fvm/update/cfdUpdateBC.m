function cfdUpdateBC(gridLevel, varargin)

if ~isempty(varargin)
    children = varargin{1};
    residual = varargin{2};
end

theCoefficients = cfdGetCoefficients;

if gridLevel==1
    theCoefficients.multigrid.BC{gridLevel,1} = theCoefficients.bc;
else    
    numberOfCoarseElements = size(children,1);
    BC = zeros(numberOfCoarseElements,1);
    
    for nCoarseElement=1:numberOfCoarseElements
        for nFineElement=1:length(children{nCoarseElement})
            iFineElement = children{nCoarseElement}(nFineElement);
            
            BC(nCoarseElement) = BC(nCoarseElement) + residual(iFineElement) ;
        end
    end    
    theCoefficients.multigrid.BC{gridLevel,1} = BC;    
end

cfdSetCoefficients(theCoefficients);