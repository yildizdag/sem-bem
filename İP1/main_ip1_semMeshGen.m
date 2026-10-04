%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Hybrid SEM-BEM
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code generates SEM-BEM mesh using the NURBS-based coarse quad
%  meshing technique.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
%
tStart = cputime;
%-Add Path:
addpath('../sem_core/')
addpath('../sem_core/geometry')
%-Read the Geometry:
FileName = 'plate_cutout_';
semPatch = [1:4]; %Enter # SEM Patches
bemPatch = []; %Enter # BEM Patches
numPatch = 4;
N = 5;
shell_dof = 6;
fluid_dof = 1;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
%
sembem2D = sembem2Dmesh(Nurbs2D,N,semPatch,bemPatch,shell_dof,fluid_dof);
%
figure('Color','w','Units','normalized','Position',[0.08 0.15 0.84 0.65]);
t = tiledlayout(1, 2,'TileSpacing','compact','Padding','compact');
%
%--------------------------------------------------------------------------
% (a) Imported NURBS geometry
%--------------------------------------------------------------------------
nexttile
%
iga2DmeshPlotNURBS(Nurbs2D);
%
axis equal
axis off
%
%--------------------------------------------------------------------------
% (b) SEM-BEM discretization
%--------------------------------------------------------------------------
nexttile
%
iga2DmeshPlotNURBS(Nurbs2D);
hold on
%
scatter3(sembem2D.nodes(:,1),sembem2D.nodes(:,2),sembem2D.nodes(:,3),60,'filled','blue');
%
hold off
%
axis equal
axis off