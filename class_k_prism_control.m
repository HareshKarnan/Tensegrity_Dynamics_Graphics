% control of the class k tensegrity prism. 
% the objective of this output feedback control is to drive the prism top
% plate to the target position
clear all;clc;close all;
% class k tensegrity dynamics with a prism example. 
[N,Cb,Cs] = tenseg_prism_sidestring(1);
B = N*Cb';S=N*Cs';b0=diag(diag(B'*B));s0=diag(diag(S'*S));
% fix the bottom 3 nodes to ground
P = [1 0 0;
     0 1 0;
     0 0 1;
     0 0 0;
     0 0 0;
     0 0 0];
D = [N(:,1) N(:,2) N(:,3)];

% external force on the structure = no external forces ! mwahaha. 
w1 = [0 0 0]';w2 = [0 0 0]';w3 = [0 0 0]';w4 = [0 0 1]';w5 = [0 0 1]';w6 = [0 0 1]';
W = 0.*[w1 w2 w3 w4 w5 w6];

% find the target node positions : 
Nt = []


