%==========================================================================
% Project             : TÜBİTAK 1001 (125M858)
% Method              : Chebyshev Spectral Element Method (CSEM)
% Meshing             : NURBS-Based Coarse-Quad Meshing
% Geometry Generation : Moving Control Point Approach
% Optimization        : Genetic Algorithm
%
% Description:
%  This is a genetic algorithm-based shape optimization framework to design 
%  stiffened plates with complex cutouts.
%
% Authors: M. Erden Yildizdag, Bekir Bediz
%==========================================================================
% BASELINE/INITIAL GEOMETRY 
clc; clear; close all;
addpath('../sem_core')
addpath('../sem_core/geometry')
% Read the Geometry imported from Rhino:
FileName1 = 'stiffOpt_test1_plate_';
numPatchPlate = 2; %Enter # Plate Patches
FileName2 = 'stiffOpt_test1_stiffener_';
numPatchStiff = 1; %Enter # Stiffener Patches
%--------------------------------------------------------------------------
% Create the Baseline/Initial Nurbs Structure (Plate+Stiffener)
%--------------------------------------------------------------------------
Nurbs2D_plate = iga2Dmesh(FileName1,numPatchPlate,1);
Nurbs2D_stiff = iga2Dmesh(FileName2,numPatchStiff,1);
%--------------------------------------------------------------------------
% Patch Connectivity (Interface: two plate patches + one stiffener patch)
%--------------------------------------------------------------------------
[Nurbs2D_plate,Nurbs2D_stiff,pconn] = patch_connectivity(Nurbs2D_plate,Nurbs2D_stiff);
%
nvars = 0;
for p = 1:size(pconn,1)
    nvars = nvars + size(Nurbs2D_plate.movingCP{pconn(p,1),pconn(p,2)},2);
end
%
dcp = 0.0.*(2.*rand(1,nvars)-1);
%
obj = evaluate_objFnc(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp);
%
figure;
hold on
iga2DmeshPlotNURBS(Nurbs2D_plate);
iga2DmeshPlotNURBS(Nurbs2D_stiff);
hold off