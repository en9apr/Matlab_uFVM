function cfdReadBCsForField(theFieldUserName)


theMesh=cfdGetMesh;
theField = cfdGetModel(theFieldUserName);
caseDirectory = theMesh.caseDirectory;
numberOfBoundaries = theMesh.numberOfBoundaries;


fieldFile = [caseDirectory,'/0/',theFieldUserName];


%
%read boundaries
%
ffid = fopen(fieldFile, 'r');
for i=1:22
    fgetl(ffid);
end
for n=1:numberOfBoundaries
        theBCUserName = fscanf(fbid, '%s', 1);
        iBC = cfdGetBoundaryIndex(theBCUserName);
        token='{';
     while(strcmp(token,'}')==false)
        token = fscanf(fbid, '%s', 1);
        if(strcmp(token,'type'))
              boundaries(iBC).type = fscanf(fbid, '%s', 1);   
        end
        if(strcmp(token,'value'))
              boundaries(iBC).startFace = fscanf(fbid, '%d', 1)+1;   
              
        end
        if(strcmp(token,'nFaces'))
              boundaries(iBC).numberOfBFaces = fscanf(fbid, '%d', 1);   
        end
        fgetl(fbid);
    end
end
mesh.boundaries = boundaries;
mesh.numberOfBoundaries = numberOfBoundaries;
mesh.numberOfPatches = numberOfBoundaries;
