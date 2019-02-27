function theMesh = processOpenFoamMesh(theMesh)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads poly mesh files from "constant/polyMesh" directory
%   and stores them in the database ~domain
%--------------------------------------------------------------------------

%
% Process basic Face Geometry
%
numberOfFaces = theMesh.numberOfFaces;
for iFace=1:numberOfFaces
    theMesh.faces(iFace).DeltaVol=0;
    iNodes = theMesh.faces(iFace).iNodes;
    numberOfiNodes = length(iNodes);
    %
    % Compute a rough centre of the face
    %
    centre = [0 0 0]';
    for iNode=iNodes
        centre = centre + theMesh.nodes(iNode).centroid;
    end
    centre = centre/numberOfiNodes;

    centroid = [0 0 0]';
    Sf = [0 0 0]';
    area = 0;
    %
    % using the centre compute the area and centoird of vitual triangles
    % based on the centre and the face nodes
    %
    for iTriangle=1:numberOfiNodes
        point1 = centre;
        point2 = theMesh.nodes(iNodes(iTriangle)).centroid;
        if(iTriangle<numberOfiNodes)
            point3 = theMesh.nodes(iNodes(iTriangle+1)).centroid;
        else 
            point3 = theMesh.nodes(iNodes(1)).centroid;
        end
        local_centroid = (point1+point2+point3)/3; 
        local_Sf  = 0.5*cross(point2-point1,point3-point1);
        local_area = cfdMagnitude(local_Sf);
 
        centroid = centroid + local_area*local_centroid;
        Sf = Sf + local_Sf;
        area = area + local_area;
    end
    centroid = centroid/area;
 
    %
    theMesh.faces(iFace).centroid = centroid;
    theMesh.faces(iFace).Sf = Sf;
    theMesh.faces(iFace).area = area;
end
%
%  compute volume and centroid of each element 
%
numberOfElements = theMesh.numberOfElements;
for iElement=1:numberOfElements
    
    iFaces = theMesh.elements(iElement).iFaces;
    %
    % Compute a rough centre of the element
    %
    centre = [0 0 0]';
    for iFace=1:length(iFaces)
        centre = centre + theMesh.faces(iFace).centroid;
    end
    centroid = [0 0 0]';
    Sf = [0 0 0]';
    centre = centre/length(iFaces);
    % using the centre compute the area and centoird of vitual triangles
    % based on the centre and the face nodes
    %
    localVolumeCentroidSum = [0 0 0]';
    localVolumeSum = 0;
    for iFace=1:length(iFaces)
        localFace = theMesh.faces(iFaces(iFace));
        localFaceSign = theMesh.elements(iElement).faceSign(iFace);
        Sf = localFace.Sf*localFaceSign;
        Cf = localFace.centroid - centre;
        % calculate face-pyramid volume
        localVolume = Sf'*Cf/3;
        % Calculate face-pyramid centre
        localCentroid = 0.75*localFace.centroid + 0.25*centre;

        %Accumulate volume-weighted face-pyramid centre
        localVolumeCentroidSum = localVolumeCentroidSum + localCentroid*localVolume;

         % Accumulate face-pyramid volume
        localVolumeSum = localVolumeSum + localVolume;
    end
    centroid = localVolumeCentroidSum/localVolumeSum; 
    volume = localVolumeSum;
    %
    theMesh.elements(iElement).volume = volume;
    theMesh.elements(iElement).OldVolume = volume;
    theMesh.elements(iElement).centroid = centroid;
end

%
% Process secondary Face Geometry
%
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
for iFace=1:numberOfInteriorFaces
    %
    theFace = theMesh.faces(iFace);
    nf = theFace.Sf/theFace.area;
    %
    element1 = theMesh.elements(theFace.iOwner);
    element2 = theMesh.elements(theFace.iNeighbour);
    %
    CN = element2.centroid - element1.centroid;
    theMesh.faces(iFace).CN = CN;
    
    eCN = CN/cfdMagnitude(CN);
    theMesh.faces(iFace).eCN = eCN;
    E = theFace.area*eCN;   
    theMesh.faces(iFace).gDiff = cfdMagnitude(E)/cfdMagnitude(CN);
    theMesh.faces(iFace).T = theFace.Sf - E;

    %
    Cf = theFace.centroid - element1.centroid;
    fF = element2.centroid - theFace.centroid;
    theMesh.faces(iFace).gf = (Cf'*nf)/(Cf'*nf + fF'*nf);
    theMesh.faces(iFace).walldist = 0;
    %
end

for iBFace=numberOfInteriorFaces+1:numberOfFaces
     %
    theBFace = theMesh.faces(iBFace);
    nf = theBFace.Sf/theBFace.area;
    %
    element1 = theMesh.elements(theBFace.iOwner);
    %
    CN = theBFace.centroid - element1.centroid;
    theMesh.faces(iBFace).CN = CN;
    theMesh.faces(iBFace).gDiff = theBFace.area* theBFace.area/dot(CN,theBFace.Sf);
    
    eCN = CN/cfdMagnitude(CN);
    theMesh.faces(iBFace).eCN = eCN;
    
    E = theBFace.area*eCN;
    
    theMesh.faces(iBFace).T = theBFace.Sf - E;
    %
    theMesh.faces(iBFace).gf = 1;
    theMesh.faces(iBFace).walldist = (CN'*theBFace.Sf)/norm(theBFace.Sf);
  
end


for iElement=1:numberOfElements
    iFaces = theMesh.elements(iElement).iFaces;
    iNeighbours = theMesh.elements(iElement).iNeighbours;
    kf=1;
    numberOfLocalInteriorFaces = length(iNeighbours);
    for iFace = iFaces(1:numberOfLocalInteriorFaces)
        if(theMesh.faces(iFace).iOwner == iElement)
            theMesh.faces(iFace).iOwnerNeighbourCoef = kf;
        elseif(theMesh.faces(iFace).iNeighbour == iElement)
            theMesh.faces(iFace).iNeighbourOwnerCoef = kf;
        end
        kf=kf+1;
    end
    
end


%loop over all the nodes and create a flag to sort boundary nodes from
%interior nodes

numberOfElements= theMesh.numberOfElements;
numberOfFaces= theMesh.numberOfFaces;
numberOfInteriorFaces= theMesh.numberOfInteriorFaces;
numberOfBoundaries=theMesh.numberOfBoundaries;
numberOfInteriorfaces=theMesh.numberOfInteriorFaces;

for iFace=1:numberOfInteriorFaces
    theMesh.faces(iFace).patchIndex=0;
 iNodes = theMesh.faces(iFace).iNodes;
 numberOfiNodes = length(iNodes);
 for iNode=1:numberOfiNodes
     Flag=1;
     theMesh.nodes(iNodes(iNode)).Flag=Flag;
 end 
end


for iBoundary=1:numberOfBoundaries
      startFace= theMesh.boundaries(iBoundary).startFace;
      numberOfBFaces=theMesh.boundaries(iBoundary).numberOfBFaces;
      
  s1=strcmp(theMesh.boundaries(iBoundary).userName,'frontAndBack');
  s2=strcmp(theMesh.boundaries(iBoundary).userName,'frontAndBackPlanes');

  if s1==1||s2==1
     for iFace=startFace:startFace+numberOfBFaces-1
 iNodes = theMesh.faces(iFace).iNodes;
 numberOfiNodes = length(iNodes);
 for iNode=1:numberOfiNodes
     Flag=1;
      theMesh.nodes(iNodes(iNode)).Flag=Flag;
 end 
     end 
  end
end


for iBoundary=1:numberOfBoundaries 
    startFace=theMesh.boundaries(iBoundary).startFace;
    numberOfBFaces=theMesh.boundaries(iBoundary).numberOfBFaces;

      s1=strcmp(theMesh.boundaries(iBoundary).userName,'frontAndBack');
  s2=strcmp(theMesh.boundaries(iBoundary).userName,'frontAndBackPlanes');

  if s1==0&&s2==0
        for iFace=startFace:startFace+numberOfBFaces-1
             theMesh.faces(iFace).patchIndex=iBoundary;
           iNodes = theMesh.faces(iFace).iNodes;
            numberOfiNodes = length(iNodes);
              for iNode=1:numberOfiNodes
                  Flag=0;
                   theMesh.nodes(iNodes(iNode)).Flag=Flag;
              end 
        end
   end
end

for iElement=1:numberOfElements
   cconn{iElement,1} = theMesh.elements(iElement).iNeighbours;
   csize(iElement,1) = length(theMesh.elements(iElement).iNeighbours);
end

theMesh.cconn =  cconn;
theMesh.csize =  csize;
end
