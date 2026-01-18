%------------------------
% SEM Axially-Loaded Bar
% Elastostatic Case
% NURBS-based Meshing
%------------------------
clc; clear; close all;
addpath('geometry')
%-Read the Geometry:
FileName = 'bar_8elm_';
numPatch = 1; %Enter #Patches
%-DOF per each control point:
local_dof = 1;
%------------------------
%-Young's Modulus
E = 200E9;
%-Geometric Props
w = 0.01;  %-width
h = 0.03;  %-height
L = 1.2;   %-length
A = w*h;   %-cross-sectional area
%-Number of Tchebychev Polynomials (per element)
N = 3;
%--------------------------------------------
% Create 1-D Nurbs Structure (reads FileName)
%--------------------------------------------
Nurbs1D = iga1Dmesh(FileName,numPatch,local_dof);
%
sem1D = sem1Dmesh(Nurbs1D,N);
% % %-Mesh
% % x = linspace(0,L,nel+1);
% % c = zeros(nel,N);
% % for i = 1:nel
% %     for j = 1:N
% %     c(i,j) = (N-1)*(i-1)+j;
% %     end
% % end
% % %-Solution
% % dof = nel*N-(nel-1);
% % K = zeros(dof,dof);
% % F = zeros(dof,1);
% % %
% % space.a=-1; space.b=1; space.N=N;
% % % calculation of necessary matrices used in the analysis
% % xi = slobat(space);
% % [F_xi,B_xi] = cheb(space);
% % [D_xi] = derivative(space);
% % [V_xi] = InnerProduct(space);
% % Q1_xi = B_xi*D_xi*F_xi;
% % %
% % xs = zeros(N,nel);
% % for i = 1:nel
% %     xs(:,i) = (0.5.*(1-xi).*x(i)+0.5.*(1+xi).*x(i+1));
% % end
% % %
% % tic;
% % for i = 1:nel
% %     c_el = c(i,:);
% %     J = (B_xi*D_xi*F_xi)*xs(:,i);
% %     Q1 = B_xi*D_xi*F_xi;
% %     kel = (E*A).*((Q1_xi./(J))'*V_xi*(Q1_xi./(J))).*(J);
% %     K(c_el,c_el) = K(c_el,c_el) + kel;
% % end
% % %
% % K(1,:) = [];
% % K(:,1) = [];
% % F(1) = [];
% % %
% % F(end) = 10000;
% % %
% % a = K\F;
% % %
% % toc;
% % a_exact = 10000*L/E/A;