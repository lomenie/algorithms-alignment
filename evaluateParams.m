dat = importdata('q1output/evaluation.tsv')

x = dat(:,1)
y = dat(:,2)
z = dat(:,5)
plot3(x,y,z,'.-')
tri = delaunay(x,y);
plot(x,y,'.')
[r,c] = size(tri);
disp(r)
h = trisurf(tri, x, y, z);
axis vis3d