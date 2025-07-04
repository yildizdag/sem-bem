% SEM-BEM LINDHOLM PLATE
% NURBS-Enhanced Coarse Quad Meshing
% Higher-Order BEM
%-------------------------------------
clc; clear; close all;
addpath('geometry')
% File name to be read:
FileName = 'sembem_rectPlate_8x8_';
semPatch = 1; %Enter # SEM Patches
bemPatch = 1; %Enter # BEM Patches
%-----------------------------------------------------------------------
% SEM mesh generator
np_u = 5; % -- FIXED!
np_v = 5; % -- FIXED!
plotNURBS = 1; % 0 or 1
plotSEM = 1; % 0 or 1
pBEM = 4;
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
%% ------------------------------------------------------------------------
% ------------------- Element Sampling ------------------------------------
% ------------------------------------------------------------------------- 
[indA,indB,indR,elementpoints,polynum_xi,polynum_eta] = element_sampling(elements,nodes);
% indA: dof(1-6)
% indR: 1 if the sampling point is shared by different elements
% elementpoints, index positions of the element points
% polynum_xi, polynum_eta: polynom numbers for elements
% -------------------------------------------------------------------------
% -------------- Assembly -----------------------------------------------
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
% ------------------------------------------------------------------------
% -------------- Eigenvalue Solution -----------------------------------
% ------------------------------------------------------------------------
tic
shift = 0.0; 
OPTS.disp = 0;
[eigVec,eigVal,flag] = eigs(Ka+shift*Ma,Ma,20,0,OPTS);
[wns,loc] = sort(real(sqrt(diag(eigVal)-shift)));
eigVec = eigVec(:,loc);
wns_HzFirst20 = (wns(1:20))/2/pi;
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
domNorm = ((sqrt(((12*(1-pois_plate^2))*(rho_plate*h_plate)/E_plate/(h_plate^3))))*10^2).*(((wns(1:20))./2./pi));
drawModeShape2; %To be fixed for multiple generation! (TUGRUL)
%
% % ind_bc1 = unique(find(abs(posn0(:,2)-max(posn0(:,2)))<1E-6));
% % ind_bc2 = unique(find(abs(posn0(:,2)-min(posn0(:,2)))<1E-6));
% % ind_bc3 = unique(find(abs(posn0(:,1)-max(posn0(:,1)))<1E-6));
% % ind_bc4 = unique(find(abs(posn0(:,1)-min(posn0(:,1)))<1E-6));
% % ind_bc = unique([ind_bc1; ind_bc2; ind_bc3; ind_bc4]);
ind_dof = transpose(setdiff(1:size(posn0,1),indAf));
U_Modes = zeros(size(posn0,1),size(UN,2));
U_Modes(ind_dof,:) = UN;
%------------------
countBEM = 0;
for i = 1:size(nodesBEM,2)
    countBEM = countBEM + size(nodesBEM{i},1);
end
H = zeros(countBEM,countBEM);
G = zeros(countBEM,countBEM);
C = 0.5.*eye(countBEM,countBEM);
modeNum=20;
b = zeros(countBEM,modeNum);
[xgp,wgp,ngp] = gaussQuad2d(12,12);
[N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),pBEM);
subd_in = [1 4 5 2
           4 7 8 5
           2 5 6 3
           5 8 9 6];
subd_out = [1 7 9 3];
%
count_col = 1;
for k=1:size(nodesBEM,2)
    %
    for j=1:size(nodesBEM{k},1)
        %
        node_i = nodesBEM{k}(j,:);
        node_ip = [node_i(1), -node_i(2), node_i(3)];
        %
        if k == 1
            ni = [0,0,1];
        elseif k == 2
            ni = [0,0,-1];
        elseif k == 3
            ni = [-1,0,0];
        elseif k == 4
            ni = [1,0,0];
        elseif k == 5
            ni = [0,-1,0];
        end
        %midpoint = elementpoints(i,13);
        %aa = find((indB==midpoint)&(indA==3));
        indd1 = find((abs(posn0(:,1)-node_i(1))<1E-6) & (abs(posn0(:,2)-node_i(2))<1E-6) & (abs(posn0(:,3)-node_i(3))<1E-6));
        if isempty(indd1)
            indd2 = find((abs(posn0(:,1)-node_i(1))<1E-6) & (abs(posn0(:,2)-node_i(2))<1E-6));
            if isempty(indd2)
                indd3 = find(abs(posn0(:,2)-node_i(2))<1E-6);
                if isempty(indd3)
                    indd4 = find((abs(posn0(:,2)-node_i(2))<1E-6));
                    if isempty(indd4)
                        indd5 = find((abs(posn0(:,1)-node_i(1))<1E-6));
                        Uix = U_Modes(indd5(1),:);
                        Uiy = U_Modes(indd5(2),:);
                        Uiz = U_Modes(indd5(3),:);
                        b(count_col,:) = ni(1).*Uix + ni(2).*Uiy + ni(3).*Uiz;
                    else
                        Uix = U_Modes(indd4(1),:);
                        Uiy = U_Modes(indd4(2),:);
                        Uiz = U_Modes(indd4(3),:);
                        b(count_col,:) = ni(1).*Uix + ni(2).*Uiy + ni(3).*Uiz;
                    end
                else
                    Uix = U_Modes(indd3(1),:);
                    Uiy = U_Modes(indd3(2),:);
                    Uiz = U_Modes(indd3(3),:);
                    b(count_col,:) = ni(1).*Uix + ni(2).*Uiy + ni(3).*Uiz;
                end
            else
                Uix = U_Modes(indd2(1),:);
                Uiy = U_Modes(indd2(2),:);
                Uiz = U_Modes(indd2(3),:);
                b(count_col,:) = ni(1).*Uix + ni(2).*Uiy + ni(3).*Uiz;
            end
        else
            Uix = U_Modes(indd1(1),:);
            Uiy = U_Modes(indd1(2),:);
            Uiz = U_Modes(indd1(3),:);
            b(count_col,:) = ni(1).*Uix + ni(2).*Uiy + ni(3).*Uiz;
        end
        %
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
                dist = norm([node_i(1)-xn(5),node_i(2)-yn(5),node_i(3)-zn(5)]);
                %
                if dist < -1
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
                    for g=1:ngp
                        %
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
%                         C(count_col,count_col) = C(count_col,count_col) + (-(1/(4*pi*r1^2))*dr_dn1)*wgp(g)*J1...
%                                                                                                                 + (-(1/(4*pi*r2^2))*dr_dn2)*wgp(g)*J2...
%                                                                                                                     + (-(1/(4*pi*r3^2))*dr_dn3)*wgp(g)*J3...
%                                                                                                                         + (-(1/(4*pi*r4^2))*dr_dn4)*wgp(g)*J4;
                        G(count_col,elementsBEM{l}(m,:)+addDOF) = G(count_col,elementsBEM{l}(m,:)+addDOF) + ((1/(4*pi*r1)-1/(4*pi*rp1))*wgp(g)*J1).*N1(1,:)...
                                                                                                                                                                                        + ((1/(4*pi*r2)-1/(4*pi*rp2))*wgp(g)*J2).*N2(1,:)...
                                                                                                                                                                                            + ((1/(4*pi*r3)-1/(4*pi*rp3))*wgp(g)*J3).*N3(1,:)...
                                                                                                                                                                                                + ((1/(4*pi*r4)-1/(4*pi*rp4))*wgp(g)*J4).*N4(1,:);
                        H(count_col,elementsBEM{l}(m,:)+addDOF) = H(count_col,elementsBEM{l}(m,:)+addDOF) + (((-1/(4*pi*r1^2))*dr_dn1+(1/(4*pi*rp1^2))*drp_dn1)*wgp(g)*J1).*N1(1,:)...
                                                                                                                                                                                        + (((-1/(4*pi*r2^2))*dr_dn2+(1/(4*pi*rp2^2))*drp_dn2)*wgp(g)*J2).*N2(1,:)...
                                                                                                                                                                                            + (((-1/(4*pi*r3^2))*dr_dn3+(1/(4*pi*rp3^2))*drp_dn3)*wgp(g)*J3).*N3(1,:)...
                                                                                                                                                                                                + (((-1/(4*pi*r4^2))*dr_dn4+(1/(4*pi*rp4^2))*drp_dn4)*wgp(g)*J4).*N4(1,:);
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
%                         C(count_col,count_col) = C(count_col,count_col) +(-(1/(4*pi*r^2))*dr_dn)*wgp(g)*J;
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
ind_bc1 = unique(find(abs(nodesBEM{1}(:,2)-max(nodesBEM{1}(:,2)))<1E-6));
ind_bc2 = unique(find(abs(nodesBEM{1}(:,2)-min(nodesBEM{1}(:,2)))<1E-6));
ind_bc3 = unique(find(abs(nodesBEM{1}(:,1)-max(nodesBEM{1}(:,1)))<1E-6));
ind_bc4 = unique(find(abs(nodesBEM{1}(:,1)-min(nodesBEM{1}(:,1)))<1E-6));
ind_bc = unique([ind_bc1; ind_bc2; ind_bc3; ind_bc4]);
C(ind_bc,ind_bc) = 0.25*eye(size(ind_bc,1),size(ind_bc,1));
%
phi = (C-H)\(G*b);
a = zeros(modeNum,modeNum);
for i = 1:modeNum
    addDOF = 0;
    for k = 1:size(elementsBEM,2)
        if k > 1
            addDOF = addDOF+size(nodesBEM{k-1},1);
        end
        for el = 1:size(elementsBEM{k},1)
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
                for j = 1:modeNum
                    a(i,j) = a(i,j) + (1000*wgp(g)*J).*(phi_g*(eigvec_g(j)));
                end
            end
        end
    end
end
A = eye(modeNum,modeNum);
c = zeros(modeNum,modeNum);
for i = 1:modeNum
    c(i,i) = eigVal(i,i);
end
[wV, wfreq] = eig(c,(a+A));
wfreq = diag(wfreq);
[wfreq2,ind] = sort((sqrt(real(wfreq))./(2*pi)));
wetV = wV(:,ind);
womNorm = ((sqrt(((12*(1-pois_plate^2))*(rho_plate*h_plate)/E_plate/(h_plate^3))))*10^2).*wfreq2;