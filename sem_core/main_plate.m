%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element Method (SEM)
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code performs free vibration analysis of isotropic Mindlin plates 
%  lie on x-y plane using a Chebyshev spectral element method combined 
%  with NURBS-based meshing for accurate geometric representation.
%
% Authors        : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'plate_20x20_';
numPatch = 1; %Enter #Patches
%-Young's Modulus
E = 200E9;
nu = 0.3;
rho = 7800;
%-Geometric Props
t = 0.01;   %-thickness
%-Number of Tchebychev Polynomials (per element)
N = 5;
modeNum = 20;
%-Element Type
ET = 1; % 1:Plate on x-y plane, 2: Plate in 3D
%-DOF per Sampling Point:
local_dof = 3;
%--------------------------------------------
% Create 1-D Nurbs Structure (reads FileName)
%--------------------------------------------
tic;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
toc;
%
tic;
sem2D = sem2Dmesh(Nurbs2D,N,local_dof);
sem2D.ET = ET;
sem2D.N = N;
sem2D.non = size(sem2D.nodes,1);
sem2D.dof = sem2D.non*local_dof;
sem2D.E = E;
sem2D.nu = nu;
sem2D.rho = rho;
sem2D.t = t;
sem2D.D = (E*t^3)/12/(1-nu^2);
sem2D.G = E/2/(1+nu);
sem2D.Ds = (5/6)*sem2D.G*t;
%
toc;
tic;
K = globalStiffness2D(sem2D);
M = globalMass2D(sem2D);
%-Boundary Conditions:
x_max = max(sem2D.nodes(:,1));
x_min = min(sem2D.nodes(:,1));
y_max = max(sem2D.nodes(:,2));
y_min = min(sem2D.nodes(:,2));
ind = find(sem2D.nodes(:,1)<x_min+1E-3|sem2D.nodes(:,1)>x_max-1E-3|sem2D.nodes(:,2)<y_min+1E-3|sem2D.nodes(:,2)>y_max-1E-3);
K([3.*ind-2;3.*ind-1;3.*ind],:) = []; K(:,[3.*ind-2;3.*ind-1;3.*ind]) = [];
M([3.*ind-2;3.*ind-1;3.*ind],:) = []; M(:,[3.*ind-2;3.*ind-1;3.*ind]) = [];
[V,freq] = eigs(K,M,modeNum,'sm');
toc;
freq = diag(sqrt(freq));
freqHZ = freq/2/pi;