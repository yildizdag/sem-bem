% SEM-BEM LINDHOLM PLATE
% NURBS-Enhanced Coarse Quad Meshing
%---------------------------------------------
clc; clear; close all;
addpath('geometry')
% File name to be read:
FileName = 'sembem_lindholmDL050_4x20_';
semPatch = 1; %Enter # SEM Patches
bemPatch = 2:6; %Enter # SEM Patches
%-----------------------------------------------------------------------
% SEM mesh generator
np_u = 5; %Sampling
np_v = 5;
plotNURBS = 1; % 0 or 1
plotSEM = 1; % 0 or 1
%-----------------------------------------------------------------------
sem2Dmesh(FileName,semPatch,bemPatch,np_u,np_v,plotNURBS,plotSEM)
%-----------------------------------------------------------------------
% DRY ANALYSIS - SEM
% clear; close all; clc
load elements
load nodes
load elementsBEM
load nodesBEM
% Material properties
E = 206.8e9;
nu = 0.3;
rhop = 7830;
h = 4.84E-3; %Thickness  0.00484             
% Boundary conditions
BCs = ['F' 'F' 'F' 'C'];    % boundary conditions, left, right, top, bottom
                                   % C: clamped, S: simply supported, 
                                   % any other letter: free
%----------------------
% Element Sampling 
%----------------------
tic
[indA,indB,indR,elementpoints,polynum_xi,polynum_eta] = element_sampling(elements,nodes);
% indA: dof(1-6)
% indR: 1 if the sampling point is shared by different elements
% elementpoints, index positions of the element points
% polynum_xi, polynum_eta: polynom numbers for elements
%-------------
% Assembly
%-------------
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
    [xelm,yelm,Kelm,Melm] = Mass_and_Stiffness_Element2(rhop,E,nu,2,2,h,polynum_xi(di1),polynum_eta(di1),locs);
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
disp(['Assembly: ' num2str(round(toc,1)) ' s'])
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
Mnorm = imag(diag(U'*(Ma)*U).^(1/2));
for i=1:size(U,2)
    UN(:,i)=U(:,i)/Mnorm(i);
end
drawModeShape2;

%ind_bc = unique([find(nodes(:,2)<0.0001); find(nodes(:,1)<0.0001); find(abs(1.4-nodes(:,2))<0.0001); find(abs(nodes(:,1)-2.0)<0.0001)]);
ind_bc = unique(find(abs(nodes(:,2)-max(nodes(:,2)))<1E-6));
ind_dof = transpose(setdiff(1:length(nodes),ind_bc));
U_Modes = zeros(5*length(nodes),size(U,2));
U_Modes([ind_dof; length(nodes)+ind_dof; 2*length(nodes)+ind_dof; 3*length(nodes)+ind_dof; 4*length(nodes)+ind_dof],:) = UN;
% %------------------
H = zeros(size(elementsBEM,1),size(elementsBEM,1));
G = zeros(size(elementsBEM,1),size(elementsBEM,1));
C = 0.5.*eye(size(elementsBEM,1));
modeNum=20;
b = zeros(size(elements,1),modeNum);
[xgp,wgp,ngp] = gaussQuad2d(8,8);
[N, dN] = linear2Dshapefun(xgp(:,1),xgp(:,2));
mid = 13;
subd_in = [1 11 13 3
                 11 21 23 13
                 13 23 25 15
                 3 13 15 5];
subd_out = [1 21 25 5];
for i=1:size(elementsBEM,1)
    disp(i)
    node_i=mean(nodesBEM(elementsBEM(i,:),:));
    node_ip = [node_i(1), -node_i(2), node_i(3)];
    a1 = nodesBEM(elementsBEM(i,3),:)-nodesBEM(elementsBEM(i,1),:);
    a2 = nodesBEM(elementsBEM(i,2),:)-nodesBEM(elementsBEM(i,1),:);
    n = cross(a1,a2)/norm(cross(a1,a2));
    %midpoint = elementpoints(i,13);
    %aa = find((indB==midpoint)&(indA==3));
    indd1 = find((abs(nodes(elements(:,13),1)-node_i(1))<1E-6) & (abs(nodes(elements(:,13),2)-node_i(2))<1E-6) & (abs(nodes(elements(:,13),3)-node_i(3))<1E-6));
    if isempty(indd1)
        indd2 = find((abs(nodes(elements(:,13),1)-node_i(1))<1E-6) & (abs(nodes(elements(:,13),2)-node_i(2))<1E-6));
        if isempty(indd2)
            indd3 = find(abs(nodes(elements(:,13),2)-node_i(2))<1E-6);
            if isempty(indd3)
                indd4 = find((abs(nodes(elements(:,11),1)-node_i(1))<1E-6));
                Uix = U_Modes(indd4(1),:);
                Uiy = U_Modes(indd4(1)+length(nodes),:);
                Uiz = U_Modes(indd4(1)+2*length(nodes),:);
                b(i,:) = n(1).*Uix + n(2).*Uiy + n(3).*Uiz;
            else
                Uix = U_Modes(indd3(1),:);
                Uiy = U_Modes(indd3(1)+length(nodes),:);
                Uiz = U_Modes(indd3(1)+2*length(nodes),:);
                b(i,:) = n(1).*Uix + n(2).*Uiy + n(3).*Uiz;
            end
        else
            Uix = U_Modes(indd2,:);
            Uiy = U_Modes(indd2+length(nodes),:);
            Uiz = U_Modes(indd2+2*length(nodes),:);
            b(i,:) = n(1).*Uix + n(2).*Uiy + n(3).*Uiz;
        end
    else
        Uix = U_Modes(indd1,:);
        Uiy = U_Modes(indd1+length(nodes),:);
        Uiz = U_Modes(indd1+2*length(nodes),:);
        b(i,:) = n(1).*Uix + n(2).*Uiy + n(3).*Uiz;
    end
    %
    % % for j=1:size(elementSBEM,1)
    % %     if j==i
    % %         for p=1:4
    % %             xn = posn0(elementpoints(j,subd_in(p,:)),1);
    % %             yn = posn0(elementpoints(j,subd_in(p,:)),2);
    % %             for k=1:ngp
    % %                 r_vector = [N(k,:)*xn; N(k,:)*yn];
    % %                 J_mat = [dN(1,:,k)*xn, dN(1,:,k)*yn; dN(2,:,k)*xn, dN(2,:,k)*yn];
    % %                 J = det(J_mat);
    % %                 r=norm(r_vector-transpose(node_i));
    % %                 rp = norm(r_vector-transpose(node_ip));
    % %                 G(i,j)=G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
    % %             end
    % %         end
    % %     else
    % %         xn = posn0(elementpoints(j,subd_out),1);
    % %         yn = posn0(elementpoints(j,subd_out),2);
    % %         for k=1:ngp
    % %             r_vector = [N(k,:)*xn; N(k,:)*yn];
    % %             J_mat = [dN(1,:,k)*xn, dN(1,:,k)*yn; dN(2,:,k)*xn, dN(2,:,k)*yn];
    % %             J = det(J_mat);
    % %             r=norm(r_vector-transpose(node_i));
    % %             rp = norm(r_vector-transpose(node_ip));
    % %             G(i,j)=G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
    % %         end
    % %     end
    % % end
end
% % phi = (C+H)\(G*b);
% % a = zeros(modeNum,modeNum);
% % for i = 1:modeNum
% %     for el = 1:size(elements,1)
% %         phi_el = phi(el,i);
% %         eigvec_el = b(el,:);
% %         for j = 1:modeNum
% %             a(i,j) = a(i,j) + (1025*phi_el*eigvec_el(j))*0.01;
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