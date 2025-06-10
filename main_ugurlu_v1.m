%-----------------------------------------------------------------------
% SEM-HOBEM RECTANGULAR VERTICAL PLATE PARTIALLY SUBMERGED INTO FLUID
% Plate in Ugurlu et al. (2008)
% NURBS-Enhaced Coarse Quad Meshing
% Higher-Order BEM -> HOBEM
%-----------------------------------------------------------------------
clc; clear; close all; format compact;
addpath('geometry')
% File name to be read:
FileName = 'sembem_rectPlate_5x5_04_';
semPatch = 1; %Enter # SEM Patches
bemPatch = 1; %Enter # BEM Patches
%-----------------------------------------------------------------------
% SEM mesh generator
np_u = 5; % -- FIXED!
np_v = 5; % -- FIXED!
plotNURBS = 1; % 0 or 1
plotSEM = 1; % 0 or 1
pBEM = 4; % 2 or 4 (quadratic or quintic)
numMode = 48;
%-----------------------------------------------------------------------
sem2Dmesh(FileName,semPatch,bemPatch,np_u,np_v,pBEM,plotNURBS,plotSEM)
%-----------------------------------------------------------------------
% DRY ANALYSIS - SEM
% clear; close all; clc
load elements
load nodes
load elementsBEM
load nodesBEM
% Material properties
rho_plate = 2400;
E_plate = 25e9;
pois_plate = 0.15;
% Dimensions, h_plate is thickness of the plate. 
h_plate = 0.15;                  % thickness of the host structure (m)
% Boundary conditions
BCs = ['C' 'C' 'C' 'C'];    % boundary conditions, left, right, top, bottom
                            % C: clamped, S: simply supported, 
                            % any other letter: free
%% ---------------------------------------------------------------------
% ------------------- Element Sampling ---------------------------------
% ----------------------------------------------------------------------
[indA,indB,indR,elementpoints,polynum_xi,polynum_eta] = element_sampling(elements,nodes);
% indA: dof(1-6)
% indR: 1 if the sampling point is shared by different elements
% elementpoints, index positions of the element points
% polynum_xi, polynum_eta: polynom numbers for elements
% -------------------------------------------------------------------------
% ---------------- Assembly -----------------------------------------------
% -------------------------------------------------------------------------
% system matrices 
Ka = sparse(size(indA,1),size(indA,1));
Ma = sparse(size(indA,1),size(indA,1));
% posn: positions for all sampling points
posn = zeros(size(indA,1),3);
for di1 = 1:size(elements,1)
    % locs: location array (25,2)
    % xlocalnow, ylocalnow: unit vector in local x- and y-directions
    [locs,xlocalnow,ylocalnow] = element_prepare1(elements(di1,:),nodes);
    % Kelm, Melm: system matrices in local domain    
    [xelm,yelm,Kelm,Melm] = Mass_and_Stiffness_Element2(rho_plate,E_plate,pois_plate,2,2,h_plate,polynum_xi(di1),polynum_eta(di1),locs);
    % indelm: indices of sampling points in the assembly matrices
    % Tnow2: transformation matrix
    [indelm,Tnow2] = element_prepare2(xlocalnow,ylocalnow,elementpoints(di1,:),indR);
    xyz0now = [0 0 0];
    posnelm = [xelm(:) yelm(:) zeros(length(xelm(:)),1)];
    zlocalnow = cross(xlocalnow,ylocalnow);    
    Tnow = [xlocalnow; ylocalnow; zlocalnow];
    posnnow = (inv(Tnow)*(posnelm'))';
    posnnow(:,1) = posnnow(:,1)+xyz0now(1);
    posnnow(:,2) = posnnow(:,2)+xyz0now(2);
    posnnow(:,3) = posnnow(:,3)+xyz0now(3);
    posnnow = round(posnnow,12);    
    posn(indelm,:) = repmat(posnnow,6,1);
    % assembly of the system matrices
    Kelmnow = [Kelm zeros(size(Kelm,1),size(Kelm,1)/5);
        zeros(size(Kelm,1)/5,size(Kelm,1)*6/5)];
    Ka(indelm,indelm) = Ka(indelm,indelm) + inv(Tnow2)*Kelmnow*Tnow2;
    Melmnow = [Melm zeros(size(Kelm,1),size(Kelm,1)/5);
        zeros(size(Kelm,1)/5,size(Kelm,1)*6/5)];
    Ma(indelm,indelm) = Ma(indelm,indelm) + inv(Tnow2)*Melmnow*Tnow2;
end
% removal of unused indices. This is mainly due to lack of one dof (rotz) 
% in fsdt modeling. 
indF = find(sum(abs(Ma))==0);
Ka(indF,:) = [];
Ka(:,indF) = [];
Ma(indF,:) = [];
Ma(:,indF) = [];
indA(indF) = [];
indB(indF) = [];
posn(indF,:) = [];
posn0 = posn;
% boundary conditions
[Ka,Ma,indA,indB,posn,indAf] = Boundary_Conditions3_plate(Ka,Ma,indA,indB,BCs,posn);
%disp(['Assembly: ' num2str(round(toc,1)) ' s'])
% ----------------------------------------------------------------------
% -------------- Eigenvalue Solution -----------------------------------
% ----------------------------------------------------------------------
tic
shift = 0.01; 
OPTS.disp = 0;
[eigVec,eigVal,flag] = eigs(Ka+shift*Ma,Ma,numMode,0,OPTS);
[wns,loc] = sort(real(sqrt(diag(eigVal)-shift)));
eigVec = eigVec(:,loc);
wns_HzFirst20 = (wns(1:numMode))/2/pi;
w1 = (wns_HzFirst20(1));
disp(['Solution time: ' num2str(round(toc,1)) ' s'])
% Sorting the mode shapes based on the sorting of the eigenvalues
U = zeros(size(eigVec,1), size(eigVec,2));
for i=1:length(eigVal)
    U(:,i) = eigVec(:,loc(i));
end
% Mass normalization of mode shapes
UN = zeros(size(U,1),size(U,2));
Mnorm = diag(U'*(Ma)*U).^(1/2);
for i=1:size(U,2)
    UN(:,i)=U(:,i)/Mnorm(i);
end
% Normalized Dry Natural Frequencies:
norm_dryFreq = ((sqrt(((12*(1-pois_plate^2))*(rho_plate*h_plate)/E_plate/(h_plate^3))))*10^2).*(((wns(1:numMode))./2./pi));
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
drawModeShape2; %To be fixed for multiple generation! (TUGRUL)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Store all displacement amplitudes:
ind_dof = transpose(setdiff(1:size(posn0,1),indAf));
U_Modes = zeros(size(posn0,1),size(UN,2));
U_Modes(ind_dof,:) = UN;
%---------------------------------------------
% Displacement Amplitudes at BEM nodes (in Z)
%---------------------------------------------
U_ModesZ_BEM = zeros(size(nodes,1),numMode);
for i = 1:numMode
    for di1 = 1:size(elements,1)
        %
        [locs,xlocalnow,ylocalnow] = element_prepare1(elements(di1,:),nodes);
        [FT_xi,IT_xi,D_xi,xi,V_xi,Q1_xi,~,space_xi] = Discretization(2, polynum_xi(di1),'xi');
        [FT_eta,IT_eta,D_eta,eta,V_eta,Q1_eta,~,space_eta] = Discretization(2, polynum_eta(di1),'eta');
        Mapping_Order = 4;
        [xelm, yelm, dxdxi, dydxi, dxdeta, dydeta, fitx, fity] = ...
            Cross_section_Mapping(Mapping_Order, locs, xi, eta);
        [indelm,Tnow2] = element_prepare2(xlocalnow,ylocalnow,elementpoints(di1,:),indR);
        deflection_w = U_Modes(indelm(51:75),i);
        FT = Fxy_mapping(polynum_xi(di1), polynum_eta(di1), FT_xi, FT_eta);
        %
        a_w = FT*deflection_w;
        %
        space_xi.N = 4;
        xi_Int = space_xi.a:(space_xi.b-space_xi.a)/space_xi.N:space_xi.b;
        space_xi.N = polynum_xi(di1);
        %
        space_eta.N = 4;
        eta_Int = space_eta.a:(space_eta.b-space_eta.a)/space_eta.N:space_eta.b;
        space_eta.N = polynum_eta(di1);
        %
        deflection_w_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, a_w);
        deflection_w_Int = reshape(transpose(reshape(deflection_w_Int,[5,5])),[25,1]);
        %
        U_ModesZ_BEM(elements(di1,:),i) = deflection_w_Int;
    end
end
%-----------------------------------------------------------------------
%------ WET ANALYSIS ---------------------------------------------------
%-----------------------------------------------------------------------
countBEM = 0;
for i = 1:size(nodesBEM,2)
    countBEM = countBEM + size(nodesBEM{i},1);
end
% Bem matrices and vectors
H = zeros(countBEM,countBEM);
G = zeros(countBEM,countBEM);
C = 0.5.*eye(countBEM,countBEM);
b = zeros(countBEM,numMode);
%------------------------------------
%------------------------------------
[xgp,wgp,ngp] = gaussQuad2d(12,12);
%------------------------------------
%------------------------------------
dist_tol = 10/20*sqrt(2);
[N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),pBEM);
count_col = 1;
for k=1:size(nodesBEM,2)
    %
    for j=1:size(nodesBEM{k},1)
        %
        node_i = nodesBEM{k}(j,:);
        node_ip = [node_i(1), -node_i(2), node_i(3)];
        %
        ni = [0,0,1];
        ind = find((abs(nodes(:,1)-node_i(1))<1E-6) & (abs(nodes(:,2)-node_i(2)))<1E-6);
        b(count_col,:) = U_ModesZ_BEM(ind,:);
        addDOF = 0;
        %
        for l=1:size(elementsBEM,2)
            %
            if l > 1
                addDOF = addDOF+size(nodesBEM{l-1},1);
            end
            for m=1:size(elementsBEM{l},1)
                %
                xn = nodesBEM{l}(elementsBEM{l}(m,:),1);
                yn = nodesBEM{l}(elementsBEM{l}(m,:),2);
                zn = nodesBEM{l}(elementsBEM{l}(m,:),3);
                %
                if pBEM == 2
                    dist = norm([node_i(1)-xn(5),node_i(2)-yn(5),node_i(3)-zn(5)]);
                elseif pBEM == 4
                    dist = norm([node_i(1)-xn(13),node_i(2)-yn(13),node_i(3)-zn(13)]);
                end
                %
                if dist < dist_tol
                    if pBEM == 2
                        sub = [2 3 6 5; 5 6 9 8; 4 5 8 7; 1 2 5 4];
                    elseif pBEM == 4
                        sub = [1 2 7 6; 2 3 8 7; 3 4 9 8; 4 5 10 9
                               6 7 12 11; 7 8 13 12; 8 9 14 13; 9 10 15 14
                               11 12 17 16; 12 13 18 17; 13 14 19 18; 14 15 20 19
                               16 17 22 21; 17 18 23 22; 18 19 24 23; 19 20 25 24];
                    end
                    if pBEM == 2
                        xi_1 = [0;0.5;1;0;0.5;1;0;0.5;1]; eta_1 = [-1;-1;-1;-0.5;-0.5;-0.5;0;0;0];
                        xi_2 = [0;0.5;1;0;0.5;1;0;0.5;1]; eta_2 = [0;0;0;0.5;0.5;0.5;1;1;1];
                        xi_3 = [-1;-0.5;0;-1;-0.5;0;-1;-0.5;0]; eta_3 = [0;0;0;0.5;0.5;0.5;1;1;1];
                        xi_4 = [-1;-0.5;0;-1;-0.5;0;-1;-0.5;0]; eta_4 = [-1;-1;-1;-0.5;-0.5;-0.5;0;0;0];
                        xn1 = [xn(2);0.5*(xn(2)+xn(3));xn(3);0.5*(xn(2)+xn(5));0.25*(xn(2)+xn(3)+xn(5)+xn(6));0.5*(xn(3)+xn(6));xn(5);0.5*(xn(5)+xn(6));xn(6)];
                        xn2 = [xn(5);0.5*(xn(5)+xn(6));xn(6);0.5*(xn(5)+xn(8));0.25*(xn(5)+xn(6)+xn(9)+xn(8));0.5*(xn(6)+xn(9));xn(8);0.5*(xn(8)+xn(9));xn(9)];
                        xn3 = [xn(4);0.5*(xn(4)+xn(5));xn(5);0.5*(xn(4)+xn(7));0.25*(xn(4)+xn(5)+xn(8)+xn(7));0.5*(xn(5)+xn(8));xn(7);0.5*(xn(7)+xn(8));xn(8)];
                        xn4 = [xn(1);0.5*(xn(1)+xn(2));xn(2);0.5*(xn(1)+xn(4));0.25*(xn(1)+xn(2)+xn(5)+xn(4));0.5*(xn(2)+xn(5));xn(4);0.5*(xn(4)+xn(5));xn(5)];
                        yn1 = [yn(2);0.5*(yn(2)+yn(3));yn(3);0.5*(yn(2)+yn(5));0.25*(yn(2)+yn(3)+yn(5)+yn(6));0.5*(yn(3)+yn(6));yn(5);0.5*(yn(5)+yn(6));yn(6)];
                        yn2 = [yn(5);0.5*(yn(5)+yn(6));yn(6);0.5*(yn(5)+yn(8));0.25*(yn(5)+yn(6)+yn(9)+yn(8));0.5*(yn(6)+yn(9));yn(8);0.5*(yn(8)+yn(9));yn(9)];
                        yn3 = [yn(4);0.5*(yn(4)+yn(5));yn(5);0.5*(yn(4)+yn(7));0.25*(yn(4)+yn(5)+yn(8)+yn(7));0.5*(yn(5)+yn(8));yn(7);0.5*(yn(7)+yn(8));yn(8)];
                        yn4 = [yn(1);0.5*(yn(1)+yn(2));yn(2);0.5*(yn(1)+yn(4));0.25*(yn(1)+yn(2)+yn(5)+yn(4));0.5*(yn(2)+yn(5));yn(4);0.5*(yn(4)+yn(5));yn(5)];
                        zn1 = [zn(2);0.5*(zn(2)+zn(3));zn(3);0.5*(zn(2)+zn(5));0.25*(zn(2)+zn(3)+zn(5)+zn(6));0.5*(zn(3)+zn(6));zn(5);0.5*(zn(5)+zn(6));zn(6)];
                        zn2 = [zn(5);0.5*(zn(5)+zn(6));zn(6);0.5*(zn(5)+zn(8));0.25*(zn(5)+zn(6)+zn(9)+zn(8));0.5*(zn(6)+zn(9));zn(8);0.5*(zn(8)+zn(9));zn(9)];
                        zn3 = [zn(4);0.5*(zn(4)+zn(5));zn(5);0.5*(zn(4)+zn(7));0.25*(zn(4)+zn(5)+zn(8)+zn(7));0.5*(zn(5)+zn(8));zn(7);0.5*(zn(7)+zn(8));zn(8)];
                        zn4 = [zn(1);0.5*(zn(1)+zn(2));zn(2);0.5*(zn(1)+zn(4));0.25*(zn(1)+zn(2)+zn(5)+zn(4));0.5*(zn(2)+zn(5));zn(4);0.5*(zn(4)+zn(5));zn(5)];
                    elseif pBEM == 4
                        subd1 = linspace(-1,-0.5,5); sub_el1x = linspace(xn(1),xn(2),5); sub_el2x = linspace(xn(2),xn(3),5); sub_el3x = linspace(xn(3),xn(4),5); sub_el4x = linspace(xn(4),xn(5),5);
                        subd2 = linspace(-0.5,0,5); sub_el1y = linspace(yn(1),yn(6),5); sub_el2y = linspace(yn(6),yn(11),5); sub_el3y = linspace(yn(11),yn(16),5); sub_el4y = linspace(yn(16),yn(21),5);
                        subd3 = linspace(0,0.5,5); sub_el1z = zeros(1,5);
                        subd4 = linspace(0.5,1,5);
                        [xi_1,eta_1] = meshgrid(subd1,subd1); xi_1 = reshape(transpose(xi_1),[25 1]); eta_1 = reshape(transpose(eta_1),[25 1]);
                        [xi_2,eta_2] = meshgrid(subd2,subd1); xi_2 = reshape(transpose(xi_2),[25 1]); eta_2 = reshape(transpose(eta_2),[25 1]);
                        [xi_3,eta_3] = meshgrid(subd3,subd1); xi_3 = reshape(transpose(xi_3),[25 1]); eta_3 = reshape(transpose(eta_3),[25 1]);
                        [xi_4,eta_4] = meshgrid(subd4,subd1); xi_4 = reshape(transpose(xi_4),[25 1]); eta_4 = reshape(transpose(eta_4),[25 1]);
                        [xi_5,eta_5] = meshgrid(subd1,subd2); xi_5 = reshape(transpose(xi_5),[25 1]); eta_5 = reshape(transpose(eta_5),[25 1]);
                        [xi_6,eta_6] = meshgrid(subd2,subd2); xi_6 = reshape(transpose(xi_6),[25 1]); eta_6 = reshape(transpose(eta_6),[25 1]);
                        [xi_7,eta_7] = meshgrid(subd3,subd2); xi_7 = reshape(transpose(xi_7),[25 1]); eta_7 = reshape(transpose(eta_7),[25 1]);
                        [xi_8,eta_8] = meshgrid(subd4,subd2); xi_8 = reshape(transpose(xi_8),[25 1]); eta_8 = reshape(transpose(eta_8),[25 1]);
                        [xi_9,eta_9] = meshgrid(subd1,subd3); xi_9 = reshape(transpose(xi_9),[25 1]); eta_9 = reshape(transpose(eta_9),[25 1]);
                        [xi_10,eta_10] = meshgrid(subd2,subd3); xi_10 = reshape(transpose(xi_10),[25 1]); eta_10 = reshape(transpose(eta_10),[25 1]);
                        [xi_11,eta_11] = meshgrid(subd3,subd3); xi_11 = reshape(transpose(xi_11),[25 1]); eta_11 = reshape(transpose(eta_11),[25 1]);
                        [xi_12,eta_12] = meshgrid(subd4,subd3); xi_12 = reshape(transpose(xi_12),[25 1]); eta_12 = reshape(transpose(eta_12),[25 1]);
                        [xi_13,eta_13] = meshgrid(subd1,subd4); xi_13 = reshape(transpose(xi_13),[25 1]); eta_13 = reshape(transpose(eta_13),[25 1]);
                        [xi_14,eta_14] = meshgrid(subd2,subd4); xi_14 = reshape(transpose(xi_14),[25 1]); eta_14 = reshape(transpose(eta_14),[25 1]);
                        [xi_15,eta_15] = meshgrid(subd3,subd4); xi_15 = reshape(transpose(xi_15),[25 1]); eta_15 = reshape(transpose(eta_15),[25 1]);
                        [xi_16,eta_16] = meshgrid(subd4,subd4); xi_16 = reshape(transpose(xi_16),[25 1]); eta_16 = reshape(transpose(eta_16),[25 1]);
                        [xn1,yn1,zn1] = meshgrid(sub_el1x,sub_el1y,sub_el1z); xn1 = reshape(transpose(xn1(:,:,1)),[25 1]); yn1 = reshape(transpose(yn1(:,:,1)),[25 1]); zn1 = reshape(transpose(zn1(:,:,1)),[25 1]);
                        [xn2,yn2,zn2] = meshgrid(sub_el2x,sub_el1y,sub_el1z); xn2 = reshape(transpose(xn2(:,:,1)),[25 1]); yn2 = reshape(transpose(yn2(:,:,1)),[25 1]); zn2 = reshape(transpose(zn2(:,:,1)),[25 1]);
                        [xn3,yn3,zn3] = meshgrid(sub_el3x,sub_el1y,sub_el1z); xn3 = reshape(transpose(xn3(:,:,1)),[25 1]); yn3 = reshape(transpose(yn3(:,:,1)),[25 1]); zn3 = reshape(transpose(zn3(:,:,1)),[25 1]);
                        [xn4,yn4,zn4] = meshgrid(sub_el4x,sub_el1y,sub_el1z); xn4 = reshape(transpose(xn4(:,:,1)),[25 1]); yn4 = reshape(transpose(yn4(:,:,1)),[25 1]); zn4 = reshape(transpose(zn4(:,:,1)),[25 1]);
                        [xn5,yn5,zn5] = meshgrid(sub_el1x,sub_el2y,sub_el1z); xn5 = reshape(transpose(xn5(:,:,1)),[25 1]); yn5 = reshape(transpose(yn5(:,:,1)),[25 1]); zn5 = reshape(transpose(zn5(:,:,1)),[25 1]);
                        [xn6,yn6,zn6] = meshgrid(sub_el2x,sub_el2y,sub_el1z); xn6 = reshape(transpose(xn6(:,:,1)),[25 1]); yn6 = reshape(transpose(yn6(:,:,1)),[25 1]); zn6 = reshape(transpose(zn6(:,:,1)),[25 1]);
                        [xn7,yn7,zn7] = meshgrid(sub_el3x,sub_el2y,sub_el1z); xn7 = reshape(transpose(xn7(:,:,1)),[25 1]); yn7 = reshape(transpose(yn7(:,:,1)),[25 1]); zn7 = reshape(transpose(zn7(:,:,1)),[25 1]);
                        [xn8,yn8,zn8] = meshgrid(sub_el4x,sub_el2y,sub_el1z); xn8 = reshape(transpose(xn8(:,:,1)),[25 1]); yn8 = reshape(transpose(yn8(:,:,1)),[25 1]); zn8 = reshape(transpose(zn8(:,:,1)),[25 1]);
                        [xn9,yn9,zn9] = meshgrid(sub_el1x,sub_el3y,sub_el1z); xn9 = reshape(transpose(xn9(:,:,1)),[25 1]); yn9 = reshape(transpose(yn9(:,:,1)),[25 1]); zn9 = reshape(transpose(zn9(:,:,1)),[25 1]);
                        [xn10,yn10,zn10] = meshgrid(sub_el2x,sub_el3y,sub_el1z); xn10 = reshape(transpose(xn10(:,:,1)),[25 1]); yn10 = reshape(transpose(yn10(:,:,1)),[25 1]); zn10 = reshape(transpose(zn10(:,:,1)),[25 1]);
                        [xn11,yn11,zn11] = meshgrid(sub_el3x,sub_el3y,sub_el1z); xn11 = reshape(transpose(xn11(:,:,1)),[25 1]); yn11 = reshape(transpose(yn11(:,:,1)),[25 1]); zn11 = reshape(transpose(zn11(:,:,1)),[25 1]);
                        [xn12,yn12,zn12] = meshgrid(sub_el4x,sub_el3y,sub_el1z); xn12 = reshape(transpose(xn12(:,:,1)),[25 1]); yn12 = reshape(transpose(yn12(:,:,1)),[25 1]); zn12 = reshape(transpose(zn12(:,:,1)),[25 1]);
                        [xn13,yn13,zn13] = meshgrid(sub_el1x,sub_el4y,sub_el1z); xn13 = reshape(transpose(xn13(:,:,1)),[25 1]); yn13 = reshape(transpose(yn13(:,:,1)),[25 1]); zn13 = reshape(transpose(zn13(:,:,1)),[25 1]);
                        [xn14,yn14,zn14] = meshgrid(sub_el2x,sub_el4y,sub_el1z); xn14 = reshape(transpose(xn14(:,:,1)),[25 1]); yn14 = reshape(transpose(yn14(:,:,1)),[25 1]); zn14 = reshape(transpose(zn14(:,:,1)),[25 1]);
                        [xn15,yn15,zn15] = meshgrid(sub_el3x,sub_el4y,sub_el1z); xn15 = reshape(transpose(xn15(:,:,1)),[25 1]); yn15 = reshape(transpose(yn15(:,:,1)),[25 1]); zn15 = reshape(transpose(zn15(:,:,1)),[25 1]);
                        [xn16,yn16,zn16] = meshgrid(sub_el4x,sub_el4y,sub_el1z); xn16 = reshape(transpose(xn16(:,:,1)),[25 1]); yn16 = reshape(transpose(yn16(:,:,1)),[25 1]); zn16 = reshape(transpose(zn16(:,:,1)),[25 1]);
                    end
                    for g=1:ngp
                        %
                        if pBEM == 4
                            ppos1 = [N(g,:)*xi_1; N(g,:)*eta_1];
                            ppos2 = [N(g,:)*xi_2; N(g,:)*eta_2];
                            ppos3 = [N(g,:)*xi_3; N(g,:)*eta_3];
                            ppos4 = [N(g,:)*xi_4; N(g,:)*eta_4];
                            ppos5 = [N(g,:)*xi_5; N(g,:)*eta_5];
                            ppos6 = [N(g,:)*xi_6; N(g,:)*eta_6];
                            ppos7 = [N(g,:)*xi_7; N(g,:)*eta_7];
                            ppos8 = [N(g,:)*xi_8; N(g,:)*eta_8];
                            ppos9 = [N(g,:)*xi_9; N(g,:)*eta_9];
                            ppos10 = [N(g,:)*xi_10; N(g,:)*eta_10];
                            ppos11 = [N(g,:)*xi_11; N(g,:)*eta_11];
                            ppos12 = [N(g,:)*xi_12; N(g,:)*eta_12];
                            ppos13 = [N(g,:)*xi_13; N(g,:)*eta_13];
                            ppos14 = [N(g,:)*xi_14; N(g,:)*eta_14];
                            ppos15 = [N(g,:)*xi_15; N(g,:)*eta_15];
                            ppos16 = [N(g,:)*xi_16; N(g,:)*eta_16];
                            [N1,~] = shapefunc2D(ppos1(1),ppos1(2),pBEM);
                            [N2,~] = shapefunc2D(ppos2(1),ppos2(2),pBEM);
                            [N3,~] = shapefunc2D(ppos3(1),ppos3(2),pBEM);
                            [N4,~] = shapefunc2D(ppos4(1),ppos4(2),pBEM);
                            [N5,~] = shapefunc2D(ppos5(1),ppos5(2),pBEM);
                            [N6,~] = shapefunc2D(ppos6(1),ppos6(2),pBEM);
                            [N7,~] = shapefunc2D(ppos7(1),ppos7(2),pBEM);
                            [N8,~] = shapefunc2D(ppos8(1),ppos8(2),pBEM);
                            [N9,~] = shapefunc2D(ppos9(1),ppos9(2),pBEM);
                            [N10,~] = shapefunc2D(ppos10(1),ppos10(2),pBEM);
                            [N11,~] = shapefunc2D(ppos11(1),ppos11(2),pBEM);
                            [N12,~] = shapefunc2D(ppos12(1),ppos12(2),pBEM);
                            [N13,~] = shapefunc2D(ppos13(1),ppos13(2),pBEM);
                            [N14,~] = shapefunc2D(ppos14(1),ppos14(2),pBEM);
                            [N15,~] = shapefunc2D(ppos15(1),ppos15(2),pBEM);
                            [N16,~] = shapefunc2D(ppos16(1),ppos16(2),pBEM);
                            posj1 = [N(g,:)*xn1; N(g,:)*yn1; N(g,:)*zn1];
                            posj2 = [N(g,:)*xn2; N(g,:)*yn2; N(g,:)*zn2];
                            posj3 = [N(g,:)*xn3; N(g,:)*yn3; N(g,:)*zn3];
                            posj4 = [N(g,:)*xn4; N(g,:)*yn4; N(g,:)*zn4];
                            posj5 = [N(g,:)*xn5; N(g,:)*yn5; N(g,:)*zn5];
                            posj6 = [N(g,:)*xn6; N(g,:)*yn6; N(g,:)*zn6];
                            posj7 = [N(g,:)*xn7; N(g,:)*yn7; N(g,:)*zn7];
                            posj8 = [N(g,:)*xn8; N(g,:)*yn8; N(g,:)*zn8];
                            posj9 = [N(g,:)*xn9; N(g,:)*yn9; N(g,:)*zn9];
                            posj10 = [N(g,:)*xn10; N(g,:)*yn10; N(g,:)*zn10];
                            posj11 = [N(g,:)*xn11; N(g,:)*yn11; N(g,:)*zn11];
                            posj12 = [N(g,:)*xn12; N(g,:)*yn12; N(g,:)*zn12];
                            posj13 = [N(g,:)*xn13; N(g,:)*yn13; N(g,:)*zn13];
                            posj14 = [N(g,:)*xn14; N(g,:)*yn14; N(g,:)*zn14];
                            posj15 = [N(g,:)*xn15; N(g,:)*yn15; N(g,:)*zn15];
                            posj16 = [N(g,:)*xn16; N(g,:)*yn16; N(g,:)*zn16];
                            a1j1 = [dN(1,:,g)*xn1; dN(1,:,g)*yn1; dN(1,:,g)*zn1];
                            a1j2 = [dN(1,:,g)*xn2; dN(1,:,g)*yn2; dN(1,:,g)*zn2];
                            a1j3 = [dN(1,:,g)*xn3; dN(1,:,g)*yn3; dN(1,:,g)*zn3];
                            a1j4 = [dN(1,:,g)*xn4; dN(1,:,g)*yn4; dN(1,:,g)*zn4];
                            a1j5 = [dN(1,:,g)*xn5; dN(1,:,g)*yn5; dN(1,:,g)*zn5];
                            a1j6 = [dN(1,:,g)*xn6; dN(1,:,g)*yn6; dN(1,:,g)*zn6];
                            a1j7 = [dN(1,:,g)*xn7; dN(1,:,g)*yn7; dN(1,:,g)*zn7];
                            a1j8 = [dN(1,:,g)*xn8; dN(1,:,g)*yn8; dN(1,:,g)*zn8];
                            a1j9 = [dN(1,:,g)*xn9; dN(1,:,g)*yn9; dN(1,:,g)*zn9];
                            a1j10 = [dN(1,:,g)*xn10; dN(1,:,g)*yn10; dN(1,:,g)*zn10];
                            a1j11 = [dN(1,:,g)*xn11; dN(1,:,g)*yn11; dN(1,:,g)*zn11];
                            a1j12 = [dN(1,:,g)*xn12; dN(1,:,g)*yn12; dN(1,:,g)*zn12];
                            a1j13 = [dN(1,:,g)*xn13; dN(1,:,g)*yn13; dN(1,:,g)*zn13];
                            a1j14 = [dN(1,:,g)*xn14; dN(1,:,g)*yn14; dN(1,:,g)*zn14];
                            a1j15 = [dN(1,:,g)*xn15; dN(1,:,g)*yn15; dN(1,:,g)*zn15];
                            a1j16 = [dN(1,:,g)*xn16; dN(1,:,g)*yn16; dN(1,:,g)*zn16];
                            a2j1 = [dN(2,:,g)*xn1; dN(2,:,g)*yn1; dN(2,:,g)*zn1];
                            a2j2 = [dN(2,:,g)*xn2; dN(2,:,g)*yn2; dN(2,:,g)*zn2];
                            a2j3 = [dN(2,:,g)*xn3; dN(2,:,g)*yn3; dN(2,:,g)*zn3];
                            a2j4 = [dN(2,:,g)*xn4; dN(2,:,g)*yn4; dN(2,:,g)*zn4];
                            a2j5 = [dN(2,:,g)*xn5; dN(2,:,g)*yn5; dN(2,:,g)*zn5];
                            a2j6 = [dN(2,:,g)*xn6; dN(2,:,g)*yn6; dN(2,:,g)*zn6];
                            a2j7 = [dN(2,:,g)*xn7; dN(2,:,g)*yn7; dN(2,:,g)*zn7];
                            a2j8 = [dN(2,:,g)*xn8; dN(2,:,g)*yn8; dN(2,:,g)*zn8];
                            a2j9 = [dN(2,:,g)*xn9; dN(2,:,g)*yn9; dN(2,:,g)*zn9];
                            a2j10 = [dN(2,:,g)*xn10; dN(2,:,g)*yn10; dN(2,:,g)*zn10];
                            a2j11 = [dN(2,:,g)*xn11; dN(2,:,g)*yn11; dN(2,:,g)*zn11];
                            a2j12 = [dN(2,:,g)*xn12; dN(2,:,g)*yn12; dN(2,:,g)*zn12];
                            a2j13 = [dN(2,:,g)*xn13; dN(2,:,g)*yn13; dN(2,:,g)*zn13];
                            a2j14 = [dN(2,:,g)*xn14; dN(2,:,g)*yn14; dN(2,:,g)*zn14];
                            a2j15 = [dN(2,:,g)*xn15; dN(2,:,g)*yn15; dN(2,:,g)*zn15];
                            a2j16 = [dN(2,:,g)*xn16; dN(2,:,g)*yn16; dN(2,:,g)*zn16];
                            J1 = norm(cross(a1j1,a2j1)); J2 = norm(cross(a1j2,a2j2)); J3 = norm(cross(a1j3,a2j3)); J4 = norm(cross(a1j4,a2j4));
                            J5 = norm(cross(a1j5,a2j5)); J6 = norm(cross(a1j6,a2j6)); J7 = norm(cross(a1j7,a2j7)); J8 = norm(cross(a1j8,a2j8));
                            J9 = norm(cross(a1j9,a2j9)); J10 = norm(cross(a1j10,a2j10)); J11 = norm(cross(a1j11,a2j11)); J12 = norm(cross(a1j12,a2j12));
                            J13 = norm(cross(a1j13,a2j13)); J14 = norm(cross(a1j14,a2j14)); J15 = norm(cross(a1j15,a2j15)); J16 = norm(cross(a1j16,a2j16));
                            nj1 = cross(a1j1,a2j1)./J1; nj2 = cross(a1j2,a2j2)./J2; nj3 = cross(a1j3,a2j3)./J3; nj4 = cross(a1j4,a2j4)./J4;
                            nj5 = cross(a1j5,a2j5)./J5; nj6 = cross(a1j6,a2j6)./J6; nj7 = cross(a1j7,a2j7)./J7; nj8 = cross(a1j8,a2j8)./J8;
                            nj9 = cross(a1j9,a2j9)./J9; nj10 = cross(a1j10,a2j10)./J10; nj11 = cross(a1j11,a2j11)./J11; nj12 = cross(a1j12,a2j12)./J12;
                            nj13 = cross(a1j13,a2j13)./J13; nj14 = cross(a1j14,a2j14)./J14; nj15 = cross(a1j15,a2j15)./J15; nj16 = cross(a1j16,a2j16)./J16;
                            r_vector1 = posj1-transpose(node_i); r_vector2 = posj2-transpose(node_i); r_vector3 = posj3-transpose(node_i); r_vector4 = posj4-transpose(node_i);
                            r_vector5 = posj5-transpose(node_i); r_vector6 = posj6-transpose(node_i); r_vector7 = posj7-transpose(node_i); r_vector8 = posj8-transpose(node_i);
                            r_vector9 = posj9-transpose(node_i); r_vector10 = posj10-transpose(node_i); r_vector11 = posj11-transpose(node_i); r_vector12 = posj12-transpose(node_i);
                            r_vector13 = posj13-transpose(node_i); r_vector14 = posj14-transpose(node_i); r_vector15 = posj15-transpose(node_i); r_vector16 = posj16-transpose(node_i);
                            r_vectorp1 = posj1-transpose(node_ip);
                            r_vectorp2 = posj2-transpose(node_ip);
                            r_vectorp3 = posj3-transpose(node_ip);
                            r_vectorp4 = posj4-transpose(node_ip);
                            r_vectorp5 = posj5-transpose(node_ip);
                            r_vectorp6 = posj6-transpose(node_ip);
                            r_vectorp7 = posj7-transpose(node_ip);
                            r_vectorp8 = posj8-transpose(node_ip);
                            r_vectorp9 = posj9-transpose(node_ip);
                            r_vectorp10 = posj10-transpose(node_ip);
                            r_vectorp11 = posj11-transpose(node_ip);
                            r_vectorp12 = posj12-transpose(node_ip);
                            r_vectorp13 = posj13-transpose(node_ip);
                            r_vectorp14 = posj14-transpose(node_ip);
                            r_vectorp15 = posj15-transpose(node_ip);
                            r_vectorp16 = posj16-transpose(node_ip);
                            r1 = norm(posj1-transpose(node_i));
                            r2 = norm(posj2-transpose(node_i));
                            r3 = norm(posj3-transpose(node_i));
                            r4 = norm(posj4-transpose(node_i));
                            r5 = norm(posj5-transpose(node_i));
                            r6 = norm(posj6-transpose(node_i));
                            r7 = norm(posj7-transpose(node_i));
                            r8 = norm(posj8-transpose(node_i));
                            r9 = norm(posj9-transpose(node_i));
                            r10 = norm(posj10-transpose(node_i));
                            r11 = norm(posj11-transpose(node_i));
                            r12 = norm(posj12-transpose(node_i));
                            r13 = norm(posj13-transpose(node_i));
                            r14 = norm(posj14-transpose(node_i));
                            r15 = norm(posj15-transpose(node_i));
                            r16 = norm(posj16-transpose(node_i));
                            rp1 = norm(posj1-transpose(node_ip));
                            rp2 = norm(posj2-transpose(node_ip));
                            rp3 = norm(posj3-transpose(node_ip));
                            rp4 = norm(posj4-transpose(node_ip));
                            rp5 = norm(posj5-transpose(node_ip));
                            rp6 = norm(posj6-transpose(node_ip));
                            rp7 = norm(posj7-transpose(node_ip));
                            rp8 = norm(posj8-transpose(node_ip));
                            rp9 = norm(posj9-transpose(node_ip));
                            rp10 = norm(posj10-transpose(node_ip));
                            rp11 = norm(posj11-transpose(node_ip));
                            rp12 = norm(posj12-transpose(node_ip));
                            rp13 = norm(posj13-transpose(node_ip));
                            rp14 = norm(posj14-transpose(node_ip));
                            rp15 = norm(posj15-transpose(node_ip));
                            rp16 = norm(posj16-transpose(node_ip));
                            dr_dn1 = (r_vector1'*nj1)/r1; dr_dn2 = (r_vector2'*nj2)/r2; dr_dn3 = (r_vector3'*nj3)/r3; dr_dn4 = (r_vector4'*nj4)/r4;
                            dr_dn5 = (r_vector5'*nj5)/r5; dr_dn6 = (r_vector6'*nj6)/r6; dr_dn7 = (r_vector7'*nj7)/r7; dr_dn8 = (r_vector8'*nj8)/r8;
                            dr_dn9 = (r_vector9'*nj9)/r9; dr_dn10 = (r_vector10'*nj10)/r10; dr_dn11 = (r_vector11'*nj11)/r11; dr_dn12 = (r_vector12'*nj12)/r12;
                            dr_dn13 = (r_vector13'*nj13)/r13; dr_dn14 = (r_vector14'*nj14)/r14; dr_dn15 = (r_vector15'*nj15)/r15; dr_dn16 = (r_vector16'*nj16)/r16;
                            drp_dn1 = (r_vectorp1'*nj1)/rp1; drp_dn2 = (r_vectorp2'*nj2)/rp2; drp_dn3 = (r_vectorp3'*nj3)/rp3; drp_dn4 = (r_vectorp4'*nj4)/rp4;
                            drp_dn5 = (r_vectorp5'*nj5)/rp5; drp_dn6 = (r_vectorp6'*nj6)/rp6; drp_dn7 = (r_vectorp7'*nj7)/rp7; drp_dn8 = (r_vectorp8'*nj8)/rp8;
                            drp_dn9 = (r_vectorp9'*nj9)/rp9; drp_dn10 = (r_vectorp10'*nj10)/rp10; drp_dn11 = (r_vectorp11'*nj11)/rp11; drp_dn12 = (r_vectorp12'*nj12)/rp12;
                            drp_dn13 = (r_vectorp13'*nj13)/rp13; drp_dn14 = (r_vectorp14'*nj14)/rp14; drp_dn15 = (r_vectorp15'*nj15)/rp15; drp_dn16 = (r_vectorp16'*nj16)/rp16;
                            G(count_col,elementsBEM{l}(m,:)+addDOF) = G(count_col,elementsBEM{l}(m,:)+addDOF) + ((1/(4*pi*r1)-1/(4*pi*rp1))*wgp(g)*J1).*N1(1,:)...
                                + ((1/(4*pi*r2)-1/(4*pi*rp2))*wgp(g)*J2).*N2(1,:)...
                                + ((1/(4*pi*r3)-1/(4*pi*rp3))*wgp(g)*J3).*N3(1,:)...
                                + ((1/(4*pi*r4)-1/(4*pi*rp4))*wgp(g)*J4).*N4(1,:)...
                                + ((1/(4*pi*r5)-1/(4*pi*rp5))*wgp(g)*J5).*N5(1,:)...
                                + ((1/(4*pi*r6)-1/(4*pi*rp6))*wgp(g)*J6).*N6(1,:)...
                                + ((1/(4*pi*r7)-1/(4*pi*rp7))*wgp(g)*J7).*N7(1,:)...
                                + ((1/(4*pi*r8)-1/(4*pi*rp8))*wgp(g)*J8).*N8(1,:)...
                                + ((1/(4*pi*r9)-1/(4*pi*rp9))*wgp(g)*J9).*N9(1,:)...
                                + ((1/(4*pi*r10)-1/(4*pi*rp10))*wgp(g)*J10).*N10(1,:)...
                                + ((1/(4*pi*r11)-1/(4*pi*rp11))*wgp(g)*J11).*N11(1,:)...
                                + ((1/(4*pi*r12)-1/(4*pi*rp12))*wgp(g)*J12).*N12(1,:)...
                                + ((1/(4*pi*r13)-1/(4*pi*rp13))*wgp(g)*J13).*N13(1,:)...
                                + ((1/(4*pi*r14)-1/(4*pi*rp14))*wgp(g)*J14).*N14(1,:)...
                                + ((1/(4*pi*r15)-1/(4*pi*rp15))*wgp(g)*J15).*N15(1,:)...
                                + ((1/(4*pi*r16)-1/(4*pi*rp16))*wgp(g)*J16).*N16(1,:);
                            H(count_col,elementsBEM{l}(m,:)+addDOF) = H(count_col,elementsBEM{l}(m,:)+addDOF) + (((-1/(4*pi*r1^2))*dr_dn1+(1/(4*pi*rp1^2))*drp_dn1)*wgp(g)*J1).*N1(1,:)...
                                + (((-1/(4*pi*r2^2))*dr_dn2+(1/(4*pi*rp2^2))*drp_dn2)*wgp(g)*J2).*N2(1,:)...
                                + (((-1/(4*pi*r3^2))*dr_dn3+(1/(4*pi*rp3^2))*drp_dn3)*wgp(g)*J3).*N3(1,:)...
                                + (((-1/(4*pi*r4^2))*dr_dn4+(1/(4*pi*rp4^2))*drp_dn4)*wgp(g)*J4).*N4(1,:)...
                                + (((-1/(4*pi*r5^2))*dr_dn5+(1/(4*pi*rp5^2))*drp_dn5)*wgp(g)*J5).*N5(1,:)...
                                + (((-1/(4*pi*r6^2))*dr_dn6+(1/(4*pi*rp6^2))*drp_dn6)*wgp(g)*J6).*N6(1,:)...
                                + (((-1/(4*pi*r7^2))*dr_dn7+(1/(4*pi*rp7^2))*drp_dn7)*wgp(g)*J7).*N7(1,:)...
                                + (((-1/(4*pi*r8^2))*dr_dn8+(1/(4*pi*rp8^2))*drp_dn8)*wgp(g)*J8).*N8(1,:)...
                                + (((-1/(4*pi*r9^2))*dr_dn9+(1/(4*pi*rp9^2))*drp_dn9)*wgp(g)*J9).*N9(1,:)...
                                + (((-1/(4*pi*r10^2))*dr_dn10+(1/(4*pi*rp10^2))*drp_dn10)*wgp(g)*J10).*N10(1,:)...
                                + (((-1/(4*pi*r11^2))*dr_dn11+(1/(4*pi*rp11^2))*drp_dn11)*wgp(g)*J11).*N11(1,:)...
                                + (((-1/(4*pi*r12^2))*dr_dn12+(1/(4*pi*rp12^2))*drp_dn12)*wgp(g)*J12).*N12(1,:)...
                                + (((-1/(4*pi*r13^2))*dr_dn13+(1/(4*pi*rp13^2))*drp_dn13)*wgp(g)*J13).*N13(1,:)...
                                + (((-1/(4*pi*r14^2))*dr_dn14+(1/(4*pi*rp14^2))*drp_dn14)*wgp(g)*J14).*N14(1,:)...
                                + (((-1/(4*pi*r15^2))*dr_dn15+(1/(4*pi*rp15^2))*drp_dn15)*wgp(g)*J15).*N15(1,:)...
                                + (((-1/(4*pi*r16^2))*dr_dn16+(1/(4*pi*rp16^2))*drp_dn16)*wgp(g)*J16).*N16(1,:);
                        elseif pBEM == 2
                            ppos1 = [N(g,:)*xi_1; N(g,:)*eta_1];
                            ppos2 = [N(g,:)*xi_2; N(g,:)*eta_2];
                            ppos3 = [N(g,:)*xi_3; N(g,:)*eta_3];
                            ppos4 = [N(g,:)*xi_4; N(g,:)*eta_4];
                            [N1,~] = shapefunc2D(ppos1(1),ppos1(2),pBEM);
                            [N2,~] = shapefunc2D(ppos2(1),ppos2(2),pBEM);
                            [N3,~] = shapefunc2D(ppos3(1),ppos3(2),pBEM);
                            [N4,~] = shapefunc2D(ppos4(1),ppos4(2),pBEM);
                            posj1 = [N(g,:)*xn1; N(g,:)*yn1; N(g,:)*zn1];
                            posj2 = [N(g,:)*xn2; N(g,:)*yn2; N(g,:)*zn2];
                            posj3 = [N(g,:)*xn3; N(g,:)*yn3; N(g,:)*zn3];
                            posj4 = [N(g,:)*xn4; N(g,:)*yn4; N(g,:)*zn4];
                            a1j1 = [dN(1,:,g)*xn1; dN(1,:,g)*yn1; dN(1,:,g)*zn1];
                            a1j2 = [dN(1,:,g)*xn2; dN(1,:,g)*yn2; dN(1,:,g)*zn2];
                            a1j3 = [dN(1,:,g)*xn3; dN(1,:,g)*yn3; dN(1,:,g)*zn3];
                            a1j4 = [dN(1,:,g)*xn4; dN(1,:,g)*yn4; dN(1,:,g)*zn4];
                            a2j1 = [dN(2,:,g)*xn1; dN(2,:,g)*yn1; dN(2,:,g)*zn1];
                            a2j2 = [dN(2,:,g)*xn2; dN(2,:,g)*yn2; dN(2,:,g)*zn2];
                            a2j3 = [dN(2,:,g)*xn3; dN(2,:,g)*yn3; dN(2,:,g)*zn3];
                            a2j4 = [dN(2,:,g)*xn4; dN(2,:,g)*yn4; dN(2,:,g)*zn4];
                            J1 = norm(cross(a1j1,a2j1));
                            J2 = norm(cross(a1j2,a2j2));
                            J3 = norm(cross(a1j3,a2j3));
                            J4 = norm(cross(a1j4,a2j4));
                            nj1 = cross(a1j1,a2j1)./J1;
                            nj2 = cross(a1j2,a2j2)./J2;
                            nj3 = cross(a1j3,a2j3)./J3;
                            nj4 = cross(a1j4,a2j4)./J4;
                            r_vector1 = posj1-transpose(node_i);
                            r_vector2 = posj2-transpose(node_i);
                            r_vector3 = posj3-transpose(node_i);
                            r_vector4 = posj4-transpose(node_i);
                            r_vectorp1 = posj1-transpose(node_ip);
                            r_vectorp2 = posj2-transpose(node_ip);
                            r_vectorp3 = posj3-transpose(node_ip);
                            r_vectorp4 = posj4-transpose(node_ip);
                            r1 = norm(posj1-transpose(node_i));
                            r2 = norm(posj2-transpose(node_i));
                            r3 = norm(posj3-transpose(node_i));
                            r4 = norm(posj4-transpose(node_i));
                            rp1 = norm(posj1-transpose(node_ip));
                            rp2 = norm(posj2-transpose(node_ip));
                            rp3 = norm(posj3-transpose(node_ip));
                            rp4 = norm(posj4-transpose(node_ip));
                            dr_dn1 = (r_vector1'*nj1)/r1;
                            dr_dn2 = (r_vector2'*nj2)/r2;
                            dr_dn3 = (r_vector3'*nj3)/r3;
                            dr_dn4 = (r_vector4'*nj4)/r4;
                            drp_dn1 = (r_vectorp1'*nj1)/rp1;
                            drp_dn2 = (r_vectorp2'*nj2)/rp2;
                            drp_dn3 = (r_vectorp3'*nj3)/rp3;
                            drp_dn4 = (r_vectorp4'*nj4)/rp4;
                            G(count_col,elementsBEM{l}(m,:)+addDOF) = G(count_col,elementsBEM{l}(m,:)+addDOF) + ((1/(4*pi*r1)-1/(4*pi*rp1))*wgp(g)*J1).*N1(1,:)...
                                + ((1/(4*pi*r2)-1/(4*pi*rp2))*wgp(g)*J2).*N2(1,:)...
                                + ((1/(4*pi*r3)-1/(4*pi*rp3))*wgp(g)*J3).*N3(1,:)...
                                + ((1/(4*pi*r4)-1/(4*pi*rp4))*wgp(g)*J4).*N4(1,:);
                            H(count_col,elementsBEM{l}(m,:)+addDOF) = H(count_col,elementsBEM{l}(m,:)+addDOF) + (((-1/(4*pi*r1^2))*dr_dn1+(1/(4*pi*rp1^2))*drp_dn1)*wgp(g)*J1).*N1(1,:)...
                                + (((-1/(4*pi*r2^2))*dr_dn2+(1/(4*pi*rp2^2))*drp_dn2)*wgp(g)*J2).*N2(1,:)...
                                + (((-1/(4*pi*r3^2))*dr_dn3+(1/(4*pi*rp3^2))*drp_dn3)*wgp(g)*J3).*N3(1,:)...
                                + (((-1/(4*pi*r4^2))*dr_dn4+(1/(4*pi*rp4^2))*drp_dn4)*wgp(g)*J4).*N4(1,:);
                        end
                    end
                else
                    for g=1:ngp
                        %
                        posj = [N(g,:)*xn; N(g,:)*yn; N(g,:)*zn];
                        a1j = [dN(1,:,g)*xn; dN(1,:,g)*yn; dN(1,:,g)*zn];
                        a2j = [dN(2,:,g)*xn; dN(2,:,g)*yn; dN(2,:,g)*zn];
                        J = norm(cross(a1j,a2j));
                        nj = cross(a1j,a2j)./J;
                        r_vector = posj-transpose(node_i);
                        r_vectorp = posj-transpose(node_ip);
                        r=norm(posj-transpose(node_i));
                        rp = norm(posj-transpose(node_ip));
                        dr_dn = (r_vector'*nj)/r;
                        drp_dn = (r_vectorp'*nj)/rp;
                        G(count_col,elementsBEM{l}(m,:)+addDOF) = G(count_col,elementsBEM{l}(m,:)+addDOF) + ((1/(4*pi*r)-1/(4*pi*rp))*wgp(g)*J).*N(g,:);
                        H(count_col,elementsBEM{l}(m,:)+addDOF) = H(count_col,elementsBEM{l}(m,:)+addDOF) + (((-1/(4*pi*r^2))*dr_dn+(1/(4*pi*rp^2))*drp_dn)*wgp(g)*J).*N(g,:);
                    end
                end
            end
        end
        count_col = count_col + 1
    end
end
%
phi = (C-H)\(G*b);
a = zeros(numMode,numMode);
for i = 1:numMode
    addDOF = 0;
    for k = 1:size(elementsBEM,2)
        %
        if k > 1
            addDOF = addDOF+size(nodesBEM{k-1},1);
        end
        for el = 1:size(elementsBEM{k},1)
            %
            phi_el = phi(elementsBEM{k}(el,:),i);
            eigvec_el = b(elementsBEM{k}(el,:),:);
            xn = nodesBEM{k}(elementsBEM{k}(el,:),1);
            yn = nodesBEM{k}(elementsBEM{k}(el,:),2);
            zn = nodesBEM{k}(elementsBEM{k}(el,:),3);
            for g=1:ngp
                %
                posj = [N(g,:)*xn; N(g,:)*yn; N(g,:)*zn];
                a1j = [dN(1,:,g)*xn; dN(1,:,g)*yn; dN(1,:,g)*zn];
                a2j = [dN(2,:,g)*xn; dN(2,:,g)*yn; dN(2,:,g)*zn];
                J = norm(cross(a1j,a2j));
                phi_g = N(g,:)*phi_el;
                eigvec_g = N(g,:)*eigvec_el;
                for j = 1:numMode
                    a(i,j) = a(i,j) + (1000*wgp(g)*J).*(phi_g*(eigvec_g(j)));
                end
            end
        end
    end
end
A = eye(numMode,numMode);
c = zeros(numMode,numMode);
for i = 1:numMode
    c(i,i) = eigVal(i,i);
end
[wV, wfreq] = eig(c,(a+A));
wfreq = diag(wfreq);
[wfreq2,ind] = sort((sqrt(real(wfreq))./(2*pi)));
wetV = wV(:,ind);
% Normalized Wet Natural Frequencies:
norm_wetFreq = ((sqrt(((12*(1-pois_plate^2))*(rho_plate*h_plate)/E_plate/(h_plate^3))))*10^2).*wfreq2;
%------------------------------------------------------------------------
%----- PLOT WET MODES ---------------------------------------------------
%------------------------------------------------------------------------
numWmode = 2;
for k = 1:numWmode
    wetModeDisp = 0.*U_Modes(:,k);
    for j = 1:numMode
        wetModeDisp = wetModeDisp - (wetV(j,k)).*U_Modes(:,j);
    end
    figure;
    hold on
    for di1 = 1:size(elements,1)
        %
        [locs,xlocalnow,ylocalnow] = element_prepare1(elements(di1,:),nodes);
        [FT_xi,IT_xi,D_xi,xi,V_xi,Q1_xi,~,space_xi] = Discretization(2, polynum_xi(di1),'xi');
        [FT_eta,IT_eta,D_eta,eta,V_eta,Q1_eta,~,space_eta] = Discretization(2, polynum_eta(di1),'eta');
        Mapping_Order = 4;
        [xelm, yelm, dxdxi, dydxi, dxdeta, dydeta, fitx, fity] = Cross_section_Mapping(Mapping_Order, locs, xi, eta);
        [indelm,Tnow2] = element_prepare2(xlocalnow,ylocalnow,elementpoints(di1,:),indR);
        deflection_w = wetModeDisp(indelm(51:75),:);
        surf(xelm,yelm,transpose(reshape(deflection_w,[5,5])));
    end
    hold off
    axis equal
    axis off
    box on
    shading interp
    colormap jet
    view(0,90)
end








