%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element Method (SEM)
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code performs free vibration analysis of Mindlin-Reissner plates 
%  using a Chebyshev spectral element method combined with NURBS-based 
%  meshing for accurate geometric representation.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'vertCompPlate_8x8_';
numPatch = 1; %Enter #Patches
%-Young's Modulus
E = 205E9;
nu = 0.3;
rho = 7800;
%-Geometric Props
t = 0.004;   %-thickness
%-Number of Tchebychev Polynomials (per element)
N = 5;
modeNum = 20;
modeNumPlot = 4;
%-Element Type:
ET = 2; % 1: Plate on x-y plane (3 DOF)
        % 2: Shell in 3D (6 DOF)
%-Formulation:
form = 2; % 1: Based on NURBS
          % 2: Based on Chebyshev
%-DOF per Sampling Point:
if ET == 1
    shell_dof = 3;
elseif ET == 2
    shell_dof = 6;
end
%----------------
% Pre-Processing
%----------------
tic;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
fprintf('NURBS data is transferred in %.4f seconds.\n', toc);
%
tic;
sem2D = sem2Dmesh(Nurbs2D,N,shell_dof);
sem2D.ET = ET;
sem2D.form = form;
sem2D.N = N;
sem2D.non = size(sem2D.nodes,1);
sem2D.dof = sem2D.non*shell_dof;
sem2D.E = E;
sem2D.nu = nu;
sem2D.rho = rho;
sem2D.t = t;
sem2D.D = (E*t^3)/12/(1-nu^2);
sem2D.G = E/2/(1+nu);
sem2D.Ds = (5/6)*sem2D.G*t;
sem2D.lame = 2*sem2D.G/(1-sem2D.nu);  
%
fprintf('SEM mesh is done in %.4f seconds.\n', toc);
tic;
%----------
% Solution
%----------
[K,M] = global2D(sem2D);
fprintf('Assembly is done in %.4f seconds.\n', toc);
%-Boundary Conditions:
tic;
%
x_min = min(sem2D.nodes(:,1)); x_max = max(sem2D.nodes(:,1));
y_min = min(sem2D.nodes(:,2)); y_max = max(sem2D.nodes(:,2));
z_min = min(sem2D.nodes(:,3));
% ind = find(sem2D.nodes(:,3)<z_min+1E-6);
ind = find(sem2D.nodes(:,1)<x_min+1E-6 | sem2D.nodes(:,1)>x_max-1E-6 |...
           sem2D.nodes(:,2)<y_min+1E-6 | sem2D.nodes(:,2)>y_max-1E-6);
BounNodes = unique([6.*ind-5; 6.*ind-4; 6.*ind-3; 6.*ind-2; 6.*ind-1; 6.*ind]);
%
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
fprintf('BCs are done in %.4f seconds.\n', toc);
tic;
%-Eigenvalue Solver
sigma = 0.5;
[V,freq] = eigs(K,M,modeNum,sigma);
[freq,loc] = sort((sqrt(diag(freq)-sigma)));
fprintf('Solution is done in %.4f seconds.\n', toc);
V = V(:,loc);
freqHz = freq/2/pi;
%
all_nodes = 1:sem2D.dof;
active = setdiff(all_nodes,BounNodes);
uModes = zeros(sem2D.dof,modeNum);
uModes(active,1:modeNum) = uModes(active,1:modeNum) + V(:,1:modeNum);
sem2D.uModes = uModes;
sem2D.freq = freq;
sem2D.freqHz = freqHz;
% % %-----------------
% % % Post-Processing
% % %-----------------
% plotModeShapes(sem2D,modeNumPlot);