%-----------------------------
% SEM Straight Timoshenko Beam
% Free Vibration
% NURBS-based Meshing
% Cantilever
%------------------------
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'bar_1000elm_';
numPatch = 1; %Enter #Patches
%-Young's Modulus
E = 200E9;
nu = 0.3;
rho = 7800;
%-Geometric Props
w = 0.5;  %-width
h = 0.03;  %-height
L = 1.2;   %-length
A = w*h;   %-cross-sectional area
%-Number of Tchebychev Polynomials (per element)
N = 3;
modeNum = 20;
%-Element Type
ET = 2; % 1:Straight Bar along x, 2: Straight Beam along x
%-DOF per Sampling Point:
local_dof = 2;
%--------------------------------------------
% Create 1-D Nurbs Structure (reads FileName)
%--------------------------------------------
Nurbs1D = iga1Dmesh(FileName,numPatch,local_dof);
%
sem1D = sem1Dmesh(Nurbs1D,N,local_dof);
sem1D.ET = ET;
sem1D.non = size(sem1D.nodes,1);
sem1D.dof = sem1D.non*local_dof;
sem1D.E = E;
sem1D.nu = nu;
sem1D.rho = rho;
sem1D.w = w;
sem1D.h = h;
sem1D.L = L;
sem1D.A = A;
sem1D.EA = E*A;
sem1D.I = h^3*w/12;
sem1D.EI = E*h^3*w/12;
sem1D.kGA = (5/6)*A*E/2/(1+nu);
%
tic;
K = globalStiffness1D(sem1D);
M = globalMass1D(sem1D);
%-Boundary Conditions:
K(1:2,:) = []; K(:,1:2) = [];
M(1:2,:) = []; M(:,1:2) = [];
[V,freq] = eigs(K,M,modeNum,'sm');
toc;
freq = diag(sqrt(freq));
freqHZ = freq/2/pi;