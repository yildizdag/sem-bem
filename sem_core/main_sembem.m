%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element - Boundary Element Method (SEM)
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code performs free vibration analysis of isotropic Mindlin plates 
%  lie on x-y plane using a Chebyshev spectral element method combined 
%  with NURBS-based meshing for accurate geometric representation.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'plate_20x20_';
numPatch = 3; %Enter #Patches
semPatch = [1 2]; %Enter # SEM Patches
bemPatch = [3]; %Enter # BEM Patches
%-Young's Modulus
E = 200E9;
nu = 0.3;
rho = 7800;
%-Geometric Props
t = 0.01;   %-thickness
%-Number of Tchebychev Polynomials (per element)
N = 5;
modeNum = 20;
modeNumPlot = 4;
%-Element Type
ET = 2; % 1:Plate on x-y plane, 2: Curved Shell
%-DOF per Sampling Point:
if ET == 1
    local_dof = 3;
elseif ET == 2
    local_dof = 6;
end
%----------------
% Pre-Processing
%----------------
tic;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
toc;
%
tic;
sembem2D = sembem2Dmesh(Nurbs2D,N,semPatch,bemPatch,shell_dof,fluid_dof);
sembem2D.ET = ET;
sembem2D.N = N;
sembem2D.non = size(sembem2D.nodes,1);
sembem2D.dof = sembem2D.non*local_dof;
sembem2D.E = E;
sembem2D.nu = nu;
sembem2D.rho = rho;
sembem2D.t = t;
sembem2D.D = (E*t^3)/12/(1-nu^2);
sembem2D.G = E/2/(1+nu);
sembem2D.Ds = (5/6)*sembem2D.G*t;
sembem2D.lame = 2*sembem2D.G/(1-sembem2D.nu);  
%
toc;
tic;
% % %----------
% % % Solution
% % %----------
% % [K,M] = global2D(sembem2D);
% % toc;
% % %-Boundary Conditions:
% % tic;
% % x_max = max(sembem2D.nodes(:,1));
% % x_min = min(sembem2D.nodes(:,1));
% % y_max = max(sembem2D.nodes(:,2));
% % y_min = min(sembem2D.nodes(:,2));
% % ind = find(sembem2D.nodes(:,1)<x_min+1E-3|sembem2D.nodes(:,1)>x_max-1E-3|sembem2D.nodes(:,2)<y_min+1E-3|sembem2D.nodes(:,2)>y_max-1E-3);
% % BounNodes = unique([6.*ind-5; 6.*ind-4; 6.*ind-3; 6.*ind-2; 6.*ind-1; 6.*ind]);
% % % BounNodes = unique([3.*ind-2; 3.*ind-1; 3.*ind]);
% % K(BounNodes,:) = []; K(:,BounNodes) = [];
% % M(BounNodes,:) = []; M(:,BounNodes) = [];
% % toc;
% % tic;
% % %-Eigenvalue Solver
% % sigma = 2000;
% % [V,freq] = eigs(K,M,modeNum,sigma);
% % [freq,loc] = sort((sqrt(diag(freq)-sigma)));
% % toc;
% % V = V(:,loc);
% % freqHz = freq/2/pi;
% % %
% % all_nodes = 1:sembem2D.dof;
% % active = setdiff(all_nodes,BounNodes);
% % uModes = zeros(sembem2D.dof,modeNum);
% % uModes(active,1:modeNum) = uModes(active,1:modeNum) + V(:,1:modeNum);
% % sembem2D.uModes = uModes;
% % sembem2D.freq = freq;
% % sembem2D.freqHz = freqHz;
% % %-----------------
% % % Post-Processing
% % %-----------------
% % % plotModeShapes(sem2D,modeNumPlot);