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
FileName = 'plate_04_';
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
toc;
%
tic;
sembem2D = sembem2Dmesh(Nurbs2D,N,semPatch,bemPatch,shell_dof,1);
sembem2D.ET = ET;
sembem2D.N = N;
sembem2D.non = size(sembem2D.nodes,1);
sembem2D.dof = sembem2D.non*shell_dof;
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
%----------
% Solution
%----------
[K,M] = global2D(sembem2D);
toc;
%-Boundary Conditions:
tic;
x_max = max(sembem2D.nodes(:,1));
x_min = min(sembem2D.nodes(:,1));
y_max = max(sembem2D.nodes(:,2));
y_min = min(sembem2D.nodes(:,2));
ind = find(sembem2D.nodes(:,1)<x_min+1E-3|sembem2D.nodes(:,1)>x_max-1E-3|sembem2D.nodes(:,2)<y_min+1E-3|sembem2D.nodes(:,2)>y_max-1E-3);
%-CCCC:
BounNodes = unique([6.*ind-5; 6.*ind-4; 6.*ind-3; 6.*ind-2; 6.*ind-1; 6.*ind]);
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
toc;
tic;
%-Eigenvalue Solver
sigma = 0;
[V,freq] = eigs(K,M,modeNum,sigma);
[freq,loc] = sort((sqrt(diag(freq)-sigma)));
toc;
V = V(:,loc);
freqHz = freq/2/pi;
%
all_nodes = 1:sembem2D.dof;
active = setdiff(all_nodes,BounNodes);
uModes = zeros(sembem2D.dof,modeNum);
uModes(active,1:modeNum) = uModes(active,1:modeNum) + V(:,1:modeNum);
sembem2D.uModes = uModes;
sembem2D.freq = freq;
sembem2D.freqHz = freqHz;
% % %-----------------
% % % Post-Processing
% % %-----------------
% % % plotModeShapes(sem2D,modeNumPlot);
% -------------------------------------
% Displacement Amplitudes at BEM nodes
% -------------------------------------
U_ModesX = zeros(size(sembem.nodes,1),modeNum);
U_ModesY = zeros(size(sembem.nodes,1),modeNum);
U_ModesZ = zeros(size(sembem.nodes,1),modeNum);
for i = 1:modeNum
    for el = 1:sembem.nel
        %
        %[locs,indelm,~] = element_prepare(sembem.elements(di1,:), sembem.nodes, sembem.elementpoints(di1,:), sembem.ind_ALL0, 0);
        %[FT_xi,IT_xi,D_xi,xi,V_xi,Q1_xi,~,space_xi] = Discretization(2, sembem.polynums(di1,1),'xi');
        %[FT_eta,IT_eta,D_eta,eta,V_eta,Q1_eta,~,space_eta] = Discretization(2, sembem.polynums(di1,2),'eta');
        %Mapping_Order = 4;
        %[xelm, yelm, dxdxi, dydxi, dxdeta, dydeta, fitx, fity] = Cross_section_Mapping(Mapping_Order, locs, xi, eta);
        el_conn = sembem2D.conn(el,:);
        deflection_u = uModes(el_conn(1:6:end),i);
        deflection_v = uModes(el_conn(2:6:end),i);
        deflection_w = uModes(el_conn(3:6:end),i);
        %FT = Fxy_mapping(sembem.polynums(di1,1), sembem.polynums(di1,2), FT_xi, FT_eta);
        %
        a_u = sembem2D.FT*deflection_u;
        a_v = sembem2D.FT*deflection_v;
        a_w = sembem2D.FT*deflection_w;
        %
        %space_xi.N = 4;
        xi_Int = linspace(-1,1,N);
        space_xi.a = -1; space_xi.b = 1;
        space_xi.N = N;
        %
        eta_Int = linspace(-1,1,N);
        space_eta.a = -1; space_eta.b = 1;
        space_eta.N = N;
        %
        deflection_u_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, a_u);
        deflection_u_Int = reshape(transpose(reshape(deflection_u_Int,[5,5])),[25,1]);
        deflection_v_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, a_v);
        deflection_v_Int = reshape(transpose(reshape(deflection_v_Int,[5,5])),[25,1]);
        deflection_w_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, a_w);
        deflection_w_Int = reshape(transpose(reshape(deflection_w_Int,[5,5])),[25,1]);
        %
        U_ModesX(sembem.elements(di1,:),i) = deflection_u_Int;
        U_ModesY(sembem.elements(di1,:),i) = deflection_v_Int;
        U_ModesZ(sembem.elements(di1,:),i) = deflection_w_Int;
    end
end
sembem.U_ModesX = U_ModesX;
sembem.U_ModesY = U_ModesY;
sembem.U_ModesZ = U_ModesZ;