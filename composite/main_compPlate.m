%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element Method (SEM)
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code performs free vibration analysis of composite Mindlin plates 
%  lie on x-y plane using a Chebyshev spectral element method combined 
%  with NURBS-based meshing for accurate geometric representation.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
%-Add Path:
addpath('../sem_core/')
addpath('geometry')
%-Read the Geometry:
FileName = 'vertCompPlate_16x16_';
numPatch = 1; %-Enter #Patches
%-Thickness:
t = 0.01;
%-Stacking Sequence:
stSeq = [0 90 90 0];
%-Material ID:
materialID = 'rectComposite';
%-Order of SEM elements:
N = 9;
%-Number of Modes to be extracted:
modeNum = 20;
%-Number of Modes to be plotted:
modeNumPlot = 4;
%-Element Type
ET = 3; % 1: Plate on x-y plane (3 DOF)
        % 2: Shell in 3D (6 DOF)
        % 3: Composite Plate on x-y plane (5 DOF)
        % 4: Composite Shell in 3D (6 DOF)
%-DOF per Sampling Point:
local_dof = 5;
%-Formulation
form = 1;
%----------------
% Pre-Processing
%----------------
tic;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
toc;
%
tic;
sem2D = sem2Dmesh(Nurbs2D,N,local_dof);
sem2D.t = t;
sem2D.stSeq = stSeq;
sem2D.ET = ET;
sem2D.local_dof = local_dof;
sem2D.form = form;
sem2D.N = N;
sem2D.non = size(sem2D.nodes,1);
sem2D.dof = sem2D.non*local_dof;
sem2D.matProp = get_material_property(materialID);
%
toc;
tic;
%----------
% Solution
%----------
[K,M] = global2D(sem2D);
toc;
%-Boundary Conditions:
tic;
x_max = max(sem2D.nodes(:,1));
x_min = min(sem2D.nodes(:,1));
y_max = max(sem2D.nodes(:,2));
y_min = min(sem2D.nodes(:,2));
%
ind1 = find(sem2D.nodes(:,1)<x_min+1E-4|sem2D.nodes(:,1)>x_max-1E-4);
ind2 = find(sem2D.nodes(:,2)<y_min+1E-4|sem2D.nodes(:,2)>y_max-1E-4);
BounNodes = unique([5.*ind1-4; 5.*ind1-3; 5.*ind1-2; 5.*ind1; 5.*ind2-4; 5.*ind2-3; 5.*ind2-2; 5.*ind2-1]);
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
toc;
tic;
%-Eigenvalue Solver
sigma = 0.0;
[V,freq] = eigs(K,M,modeNum,'sm');
[freq,loc] = sort((sqrt(diag(freq)-sigma)));
toc;
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
%-----------------
% Post-Processing
%-----------------
%plotModeShapes(sem2D,modeNumPlot);