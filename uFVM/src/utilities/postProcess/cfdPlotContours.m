function cfdPlotContours(phi, varargin)


opt = struct('N',    20, ...
                'min', [], ...
                'max', [], ...
                'mesh', [], ...
                'cmap',  @jet,...
                'basealpha', 1, ...
                'binc',    [],...
                'patchn', 1000000,...
                'interpolant', []);
            
            
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfFaces = theMesh.numberOfFaces;
gec = [theMesh.elements.centroid];
gfc = [theMesh.faces(numberOfInteriorFaces+1:numberOfFaces).centroid];
gc=[gec gfc]';

Mgc = max(gc);
mgc = min(gc);

cmesh = [20  20 20];
tmp = cell(3, 1);
for i = 1:3
    m = mgc(i);
    M = Mgc(i);
    tmp{i} = m:((M-m)/cmesh(i)):M;
end


[x,y,z] = meshgrid(tmp{:});
%F = scatteredInterpolant(gc(:,1), gc(:,2), gc(:,3), phi);
F=@(xi, yi, zi)griddata(gc(:,1), gc(:,2), gc(:,3), phi, xi, yi, zi);
interpolant = F(x,y,z);
N=20;
binc = min(phi):(max(phi)-min(phi))/(N-1):max(phi);
c = opt.cmap(N);
figure(4)
clf;
hold on;

%axis off
axis equal
for i = 1:N
   patch(isosurface(x, y, z, interpolant, binc(i)), 'facec', c(i,:), 'edgec', 'none', 'facea', opt.basealpha*i/(2*N));
end
cfdPlotMesh('FaceAlpha', 0, 'EdgeAlpha', .1)

end
