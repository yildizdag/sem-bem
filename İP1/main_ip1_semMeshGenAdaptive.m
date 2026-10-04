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
addpath('geometry')
%-Read the Geometry:
FileName = 'ip1_stiffenedPlate_';
semPatch = [1:3]; %Enter # SEM Patches
bemPatch = []; %Enter # BEM Patches
numPatch = length(semPatch)+length(bemPatch);
%
shell_dof = 6;
fluid_dof = 1;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
%
delta0 = 0.15;
ls = 1;
sembem2D = sembem2DmeshAdaptive(Nurbs2D,delta0,ls,semPatch,bemPatch,shell_dof,fluid_dof);
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
view(-35,25)
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
view(-35,25)