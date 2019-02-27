function cfdSetFixedElement(varargin)

global Domain;

if nargin==0
    algorithm = cfdGetAlgorithm;
    iElement = Domain.foam.fvSolution.(algorithm).pRefCell;
else
    iElement = varargin{1};
end

Domain.iFixedElement = iElement;
