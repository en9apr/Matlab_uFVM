function dphi = cfdSolveILU(ac,anb,bc,dc,rc,cconn,dphi)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

% ILU Iterate
%
%
numberOfElements = length(ac);
%..........................
% update Residual array
%..........................
for iElement=1:numberOfElements
    % Compute Residual
    conn = cconn{iElement};
    res = -ac(iElement)*dphi(iElement) + bc(iElement);
    theNumberOfNeighbours = length(conn);
    %
    for iLocalNeighbour=1:theNumberOfNeighbours
        % get the neighbour cell index
        j = conn(iLocalNeighbour);
        res = res - anb{iElement}(iLocalNeighbour)*dphi(j);
    end
    %
    rc(iElement) = res;
end
%..........................
% Forward Substitution
%..........................
for i1=1:numberOfElements
    %
    %
    mat1 = dc(i1)*rc(i1);
    %
    i1NNb = length(cconn{i1});
    i1NbList = cconn{i1};
    
    j1_ = 0;
    % loop over neighbours of i
    while((j1_+1) <= i1NNb)
        j1_ = j1_ +1;
        j1 = i1NbList(j1_);
        % for all neighbour j > i do
        if(j1 > i1 && j1<=numberOfElements)
            %if(j>i){    //*** The old sequential version, left for backward debugging
            j1NbList = cconn{j1};
            j1NNB = length(j1NbList);
            
            i1_= 0;
            k = 0;
            % get A[j][i]
            while( ((i1_+1)<=j1NNB) && (k ~= i1))
                i1_ = i1_ +1;
                k = j1NbList(i1_);
            end
            % Compute rc
            if(k == i1 )
                mat2 =  anb{j1}(i1_)*mat1;
                rc(j1) = rc(j1) - mat2;
            else
                disp('ILU Solver Error The index for i  in element j  is not found \n');
            end
        end
    end
end

% backward substitution
for i1=numberOfElements:-1:1
    %
    if(i1 < (numberOfElements))
        i1NBList = cconn{i1};
        i1NNb = length(i1NBList);
        j1_ = 0;
        % loop over neighbours of i
        while((j1_+1) <= i1NNb)
            j1_ = j1_ +1;
            j = i1NBList(j1_);
            if(j>i1)
                rc(i1) = rc(i1) - anb{i1}(j1_)*rc(j);
            end
        end
    end
    % compute product D[i]*R[i]
    mat1 = dc(i1)*rc(i1);
    rc(i1) = mat1;
    % update dphi
    dphi(i1) = dphi(i1) + mat1;
end

