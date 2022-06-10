
clear all;
close all;

u = load('u.mat','-ascii');
tri = load('original_mesh_tri.mat','-ascii');
x = load('original_mesh_x.mat','-ascii');
y = load('original_mesh_y.mat','-ascii');

trisurf(tri,x,y,u);