function fvmNodes=cfdSetupNodePatches(fvmFaces,numberOfNodes,fvmNodes)


exteriorNodes=[];
interiorNodes=[];


exteriorFacesIndices = find([fvmFaces(:).patchIndex]~=0);

for iBExteriorFaces=1:length(exteriorFacesIndices)
    exteriorNodes=[exteriorNodes fvmFaces(exteriorFacesIndices(iBExteriorFaces)).nodes];
end

interiorFacesIndices = find([fvmFaces(:).patchIndex]==0);

for iInteriorFaces=1:length(interiorFacesIndices)
    interiorNodes=[interiorNodes fvmFaces(interiorFacesIndices(iInteriorFaces)).nodes];
end

%remove redundant nodes
red=[];
red1=[];
interiorNodesfin=[];
Flag=0;

for i=2:length(exteriorNodes)
    for j=1:i-1
        if (exteriorNodes(i)==exteriorNodes(j))
            red1=[red1 i];
            break
        end
    end
end
for i=1:length(red1)
    exteriorNodes(red1(i))=[];
    for j=i:length(red1)
        red1(j)=red1(j)-1;
    end
end


% for i=2:length(interiorNodes)
%     for j=1:i-1
%         if (interiorNodes(i)==interiorNodes(j))
%             red=[red i];
%             break
%         end
%     end
% end
% for i=1:length(red)
%     interiorNodes(red(i))=[];
%     for j=i:length(red)
%         red(j)=red(j)-1;
%     end
% end
% 
% for i=1:length(interiorNodes)
%     Flag=0;
%     for j=1:length(exteriorNodes)
%         if(interiorNodes(i)==exteriorNodes(j))
%             Flag=Flag+1;
%             break
%           
%         end
%     end
%     if (Flag==0)
%        interiorNodesfin=[interiorNodesfin interiorNodes(i)];
%     end
% end

for i=1:numberOfNodes
    Flag2=0;
    for j=1:length(exteriorNodes)
        if(fvmNodes(i).id==exteriorNodes(j))
            Flag2=1;
            
            break
            
        end
    end
    if(Flag2~=0)
        fvmNodes(i).patchIndex=1;
    else
        fvmNodes(i).patchIndex=0;
    end
end

        
            