% SEM-BEM LINDHOLM PLATE
% NURBS-Enhanced Coarse Quad Meshing
% Higher-Order BEM
%-------------------------------------
clc; clear; close all;
addpath('geometry')
% File name to be read:
FileName = 'sembem_lindholmDL050_4x20_';
semPatch = 1; %Enter # SEM Patches
bemPatch = 2:6; %Enter # BEM Patches
%-----------------------------------------------------------------------
% SEM mesh generator
np_u = 5; % -- FIXED!
np_v = 5; % -- FIXED!
plotNURBS = 1; % 0 or 1
plotSEM = 1; % 0 or 1
pBEM = 2;
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
rho_plate = 7830;
E_plate = 206.8e9;
pois_plate = 0.3;
% Dimensions, h_plate is thickness of the plate. 
h_plate = 0.00484;                  % thickness of the host structure (m)
% Boundary conditions
BCs = ['F' 'F' 'F' 'C'];    % boundary conditions, left, right, top, bottom
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
[Ka,Ma,indA,indB,posn] = Boundary_Conditions3_plate(Ka,Ma,indA,indB,BCs,posn);
%disp(['Assembly: ' num2str(round(toc,1)) ' s'])
% ------------------------------------------------------------------------
% -------------- Eigenvalue Solution -----------------------------------
% ------------------------------------------------------------------------
tic
shift = 0.0; 
OPTS.disp=0;
[eigVec,eigVal,flag] = eigs(Ka+shift*Ma,Ma,20,0,OPTS);
[wns,loc] = sort(real(sqrt(diag(eigVal)-shift)));
eigVec = eigVec(:,loc);
wns_HzFirst10 = (wns(1:20))/2/pi;
w1 = (wns_HzFirst10(1));
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
drawModeShape2; %To be fixed for multiple generation! (TUGRUL)
%
% % ind_bc = unique(find(abs(posn0(:,2)-max(posn0(:,2)))<1E-6));
% % ind_dof = transpose(setdiff(1:size(posn0,1),ind_bc));
% % U_Modes = zeros(size(posn0,1),size(UN,2));
% % U_Modes(ind_dof,:) = UN;
% % %------------------
% % countBEM = 0;
% % for i = 1:size(nodesBEM,2)
% %     countBEM = countBEM + size(nodesBEM{i},1);
% % end
% % H = zeros(countBEM,countBEM);
% % G = zeros(countBEM,countBEM);
% % C = 0.5.*eye(countBEM);
% % modeNum=20;
% % b = zeros(countBEM,modeNum);
% % [xgp,wgp,ngp] = gaussQuad2d(8,8);
% % [N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),pBEM);
% % subd_in = [1 4 5 2
% %                  4 7 8 5
% %                  2 5 6 3
% %                  5 8 9 6];
% % subd_out = [1 7 9 3];
% % cel = [1 8 4 5 9 7 2 6 3];
% % count_col = 1;
% % for k=1:size(nodesBEM,2)
% %     %
% %     for j=1:size(nodesBEM{k},1)
% %         %
% %         node_i = nodesBEM{k}(j,:);
% %         node_ip = [node_i(1), -node_i(2), node_i(3)];
% %         %
% %         addDOF = 0;
% %         for l=1:size(elementsBEM,2)
% %             %
% %             for m=1:size(elementsBEM{l},1)
% %                 %
% %                 xn = nodesBEM{l}(elementsBEM{l}(m,cel),1);
% %                 yn = nodesBEM{l}(elementsBEM{l}(m,cel),2);
% %                 zn = nodesBEM{l}(elementsBEM{l}(m,cel),3);
% %                 for g=1:ngp
% %                     %
% %                     posj = [N(g,:)*xn; N(g,:)*yn; N(g,:)*zn];
% %                     a1j = [dN(1,:,g)*xn; dN(1,:,g)*yn; dN(1,:,g)*zn];
% %                     a2j = [dN(2,:,g)*xn; dN(2,:,g)*yn; dN(2,:,g)*zn];
% %                     J = norm(cross(a1j,a2j));
% %                     nj = cross(a1j,a2j)./J;
% %                     r_vector = posj-transpose(node_i);
% %                     r_vectorp = posj-transpose(node_ip);
% %                     r=norm(posj-transpose(node_i));
% %                     rp = norm(posj-transpose(node_ip));
% %                     dr_dn = (r_vector'*nj)/r;
% %                     drp_dn = (r_vectorp'*nj)/rp;
% %                     G(count_col,elementsBEM{l}(m,cel)+addDOF) = G(count_col,elementsBEM{l}(m,cel)+addDOF) + ((1/(4*pi*r)-1/(4*pi*rp))*wgp(g)*J).*N(g,:);
% %                     H(count_col,elementsBEM{l}(m,cel)+addDOF) = H(count_col,elementsBEM{l}(m,cel)+addDOF) + (((-1/(4*pi*r^2))*dr_dn+(1/(4*pi*rp^2))*drp_dn)*wgp(g)*J).*N(g,:);
% %                 end
% %             end
% %             if l > 1
% %                 addDOF = addDOF+size(nodesBEM{l-1},1);
% %             end
% %         end
% %         count_col = count_col + 1;
% %     end
% % end


% %         if j==i
% %             for p=1:4
% %                 xn = nodesBEM(elementsBEM(j,subd_in(p,:)),1);
% %                 yn = nodesBEM(elementsBEM(j,subd_in(p,:)),2);
% %                 zn = nodesBEM(elementsBEM(j,subd_in(p,:)),3);
% %                 for k=1:ngp
% %                     posj = [N(k,:)*xn; N(k,:)*yn; N(k,:)*zn];
% %                     a1j = [dN(1,:,k)*xn; dN(1,:,k)*yn; dN(1,:,k)*zn];
% %                     a2j = [dN(2,:,k)*xn; dN(2,:,k)*yn; dN(2,:,k)*zn];
% %                     J = norm(cross(a1j,a2j));
% %                     r=norm(posj-transpose(node_i));
% %                     rp = norm(posj-transpose(node_ip));
% %                     G(i,j)=G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
% %                     H(i,j) = H(i,j)+((-1/(4*pi*r^2))*dr_dn+(1/(4*pi*rp^2))*drp_dn)*wgp(g)*J;
% %                 end
% %             end
% %         else
% %             xn = nodesBEM(elementsBEM(j,subd_out),1);
% %             yn = nodesBEM(elementsBEM(j,subd_out),2);
% %             zn = nodesBEM(elementsBEM(j,subd_out),3);
% %             for k=1:ngp
% %                 posj = [N(k,:)*xn; N(k,:)*yn; N(k,:)*zn];
% %                 a1j = [dN(1,:,k)*xn; dN(1,:,k)*yn; dN(1,:,k)*zn];
% %                 a2j = [dN(2,:,k)*xn; dN(2,:,k)*yn; dN(2,:,k)*zn];
% %                 J = norm(cross(a1j,a2j));
% %                 nj = cross(a1j,a2j)./J;
% %                 r_vector = posj-transpose(node_i);
% %                 r_vectorp = posj-transpose(node_ip);
% %                 r=norm(r_vector);
% %                 rp = norm(r_vectorp);
% %                 dr_dn = (r_vector'*nj)/r;
% %                 drp_dn = (r_vectorp'*nj)/rp;
% %                 G(i,j) = G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
% %                 H(i,j) = H(i,j)+((-1/(4*pi*r^2))*dr_dn+(1/(4*pi*rp^2))*drp_dn)*wgp(k)*J;
% %             end
% %         end
% %     end
% % end
% % phi = (C-H)\(G*b);
% % a = zeros(modeNum,modeNum);
% % for i = 1:modeNum
% %     for el = 1:size(elementsBEM,1)
% %         phi_el = phi(el,i);
% %         eigvec_el = b(el,:);
% %         l1_el = norm(nodesBEM(elementsBEM(el,7),:)-nodesBEM(elementsBEM(el,1),:));
% %         l2_el = norm(nodesBEM(elementsBEM(el,3),:)-nodesBEM(elementsBEM(el,1),:));
% %         A_el = l1_el*l2_el;
% %         for j = 1:modeNum
% %             a(i,j) = a(i,j) + A_el*(1025*phi_el*eigvec_el(j));
% %         end
% %     end
% % end
% % A = eye(modeNum,modeNum);
% % c = zeros(modeNum,modeNum);
% % for i = 1:modeNum
% %     c(i,i) = eigVal(i,i);
% % end
% % [wV, wfreq] = eig(c,(a+A));
% % wfreq = diag(wfreq);
% % [wfreq2,ind] = sort((sqrt(real(wfreq))./(2*pi)));
% % wetV = wV(:,ind);