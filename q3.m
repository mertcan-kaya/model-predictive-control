clc, clear all, close all

% State-pace parameters
A = 0.9;
B = 0.5;
C = 1;
D = 0;

Q = 1;
R = 0;

Np = 30;
Nu = 2;

%% a)

A_aug = [A B;0 1];
B_aug = [B;0];
C_aug = [C 0];

%% b)

Lx = 0;
Ly = 0;

L = [Lx;Ly];



