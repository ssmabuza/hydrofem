clear all;
close all;
clc;

u = load('conc.mat','--ascii');
x = load('original_mesh_x.mat','--ascii');
y = load('original_mesh_y.mat','--ascii');
tri = load('original_mesh_tri.mat','--ascii');

time_ind = 3;
delta_t = 0.03125;
time = time_ind*delta_t;

% figure from hydrofem
fig1 = figure(1)
trisurf(tri,x,y,u(time_ind,:));
view(2);
colormap(jet);
colorbar();
axis tight;
axis equal;


