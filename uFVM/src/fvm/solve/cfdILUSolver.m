function output = cfdILUSolver(option)


%===================================================
%  written by the CFD Group @ AUB, Fall 2014
%===================================================
if(strcmp(option,'setup'))
    output = cfdILUFactorize;
else
    output = cfdILUIterate(10);
end
end

%===================================================
%  cfdILUFactorize
%===================================================

function theCoefficients =  cfdILUFactorize
%
theCoefficients = cfdGetCoefficients;

ac  = theCoefficients.ac;
anb = theCoefficients.anb;
bc  = theCoefficients.bc;
numberOfElements = length(ac);

dc = theCoefficients.dc;
rc = theCoefficients.rc;

for i1=1:numberOfElements
    dc(i1) = ac(i1);
end

for i1=1:numberOfElements
    dc(i1) = 1.0/dc(i1);
    rc(i1) = bc(i1);
    
    i1NbList = theCoefficients.cconn{i1};
    i1NNb = length(i1NbList);
    
    if(i1~=numberOfElements-1)
        % loop over neighbours of iElement
        j1_ = 1;
        while(j1_<=i1NNb)            
            jj1 = i1NbList(j1_);
            % for all neighbour j > i do
            if((jj1>i1) && (jj1<=numberOfElements))
                j1NbList = theCoefficients.cconn{jj1};
                j1NNb = length(j1NbList);
                i1_= 0;
                k1 = -1;
                % find _i index to get A[j][_i]
                while((i1_<=j1NNb) && (k1 ~= i1))
                    i1_ = i1_ + 1;
                    k1 = j1NbList(i1_);
                end
                % Compute A[j][i]*D[i]*A[i][j]
                if(k1 == i1)
                    dc(jj1) = dc(jj1) - anb{jj1}(i1_)*dc(i1)*anb{i1}(j1_);
                else
                    disp('the index for i in j is not found');
                end
            end
            j1_ = j1_ + 1;
        end
    end
end
theCoefficients.dc = dc;
theCoefficients.rc = rc;
theCoefficients.dphi = 0*theCoefficients.dphi;
cfdSetCoefficients(theCoefficients);

end

%===================================================
%  cfdILUIterate
%===================================================


function Residual0 = cfdILUIterate(theNumberOfIterations)

theCoefficients = cfdGetCoefficients;

ac = theCoefficients.ac;
anb = theCoefficients.anb;
bc = theCoefficients.bc;
rc = theCoefficients.rc;
dc = theCoefficients.dc;
dphi = theCoefficients.dphi;
%
theNumberOfElements = length(bc);
%
%
for NIter=1:theNumberOfIterations;
    %..........................
    % update Residual array
    %..........................
    for iElement=1:theNumberOfElements
        % Compute Residual
        conn = theCoefficients.cconn{iElement};
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
    %
    Residual1 = sqrt(sum(rc.*rc)/theNumberOfElements);
    if(NIter==1)
        Residual0 = Residual1;
    end
    %..........................
    % Forward Substitution
    %..........................
    for i1=1:theNumberOfElements
        %
        %
        mat1 = dc(i1)*rc(i1);
        %
        i1NNb = length(theCoefficients.cconn{i1});
        i1NbList = theCoefficients.cconn{i1};
        
        j1_ = 0;
        % loop over neighbours of i
        while((j1_+1) <= i1NNb)
            j1_ = j1_ +1;
            j1 = i1NbList(j1_);
            % for all neighbour j > i do
            if(j1 > i1 && j1<=theNumberOfElements)
                %if(j>i){    //*** The old sequential version, left for backward debugging
                j1NbList = theCoefficients.cconn{j1};
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
    %
    % .........................
    % backward substitution
    % ..........................
    for i1=theNumberOfElements:-1:1
        %
        if(i1 < (theNumberOfElements))
            i1NBList = theCoefficients.cconn{i1};
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
    %         theCoefficients.rc = rc;
    %         theCoefficients.dphi = dphi;
    %         cfdSetCoefficients(theCoefficients);
    %         theCoefficients = cfdGetCoefficients;
end
theCoefficients.rc = rc;
theCoefficients.dphi = dphi;
cfdSetCoefficients(theCoefficients);
end
%
%
