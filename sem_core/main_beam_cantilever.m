%-----------------------------
% SEM Straight Timoshenko Beam
% Free Vibration
% NURBS-based Meshing
% Cantilever
%------------------------
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'bar_8elm_';
numPatch = 1; %Enter #Patches
%-Young's Modulus
E = 200E9;
%-Geometric Props
w = 0.01;  %-width
h = 0.03;  %-height
L = 1.2;   %-length
A = w*h;   %-cross-sectional area
%-Number of Tchebychev Polynomials (per element)
N = 5;
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
sem1D.w = w;
sem1D.h = h;
sem1D.L = L;
sem1D.A = A;
sem1D.EA = E*A;
%
tic;
K = globalStiffness1D(sem1D);
F = zeros(sem1D.dof,1);
%
K(1,:) = [];
K(:,1) = [];
F(1) = [];
%
F(end) = 10000;
%
a = K\F;
%
toc;
a_exact = 10000*L/E/A;