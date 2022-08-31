clear all;
close all;
clc;

u = load('u.mat','--ascii');
x = load('original_mesh_x.mat','--ascii');
y = load('original_mesh_y.mat','--ascii');
tri = load('original_mesh_tri.mat','--ascii');

% figure from hydrofem
figure(1)
trisurf(tri,x,y,u);

% plot the exact solution
xx = 0:0.02:1.0;
yy = xx;
[X, Y] = meshgrid(xx,yy);
T = delaunay(X,Y);
U = sin(2.0*pi*X).*sin(2.0*pi*Y);

% exact solution figure
figure(2)
trisurf(T,X,Y,U)
