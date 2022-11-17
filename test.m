rng('default')
xy = -2.5 + 5*rand([200 2]);
x = xy(:,1);
y = xy(:,2);
v = x.*exp(-x.^2-y.^2);

[x, y, z] = ptCloud2Mesh(x, y, v, 0.01, 0.01)