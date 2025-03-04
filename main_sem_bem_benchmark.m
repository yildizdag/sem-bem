% NURBS-Enhanced Mesh Generator for
% Spectral Element Method (SEM)
%---------------------------------------
clc; clear; close all;
addpath('geometry')
% Read the Geometry imported from Rhino:
FileName = 'sembem_lindholmDL050_4x20_';
numPatch = 1; %Enter # Patches
% Degrees of Freedom per each node:
local_dof = 1;
%---------------------------------------
% Create 2-D Nurbs Structure (reads FileName)
%---------------------------------------
Nurbs2D = iga2Dmesh(FileName,numPatch,local_dof);
%---------------------------------------
% Refinement (if necessary)
%---------------------------------------
% ur = 0; % Refinement Level in u direction
% vr = 0; % Refinement Level in v direction
% Nurbs2D = hrefine2D(Nurbs2D,1,ur,vr);
% Nurbs2D = hrefine2D(Nurbs2D,2,ur,vr);
% Nurbs2D = hrefine2D(Nurbs2D,3,ur,vr);
% Nurbs2D = hrefine2D(Nurbs2D,4,ur,vr);
% Nurbs2D = hrefine2D(Nurbs2D,5,ur,vr);
%---------------------------------------
% Plot Imported 2-D NURBS Structure
%---------------------------------------
figure;
iga2DmeshPlotNURBS(Nurbs2D);
axis off
%-------------------------------------------------
% Points for Spectral Element Method (e.g. 5 x 5, 3 x 3, etc.)
%-------------------------------------------------
np_u = 5;
np_v = 5;
tot_el = 0; %Total num of elements
for k = 1:Nurbs2D.numpatch
    tot_el = tot_el + Nurbs2D.nel{k};
end
elData = zeros(np_u*np_v,3,tot_el);
nodeData = zeros(np_u*np_v*tot_el,3);
count_el = 1;
count_node = 1;
figure;
iga2DmeshPlotNURBS(Nurbs2D);
for k = 1:Nurbs2D.numpatch
    for el = 1:Nurbs2D.nel{k}
        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
        u_sample = linspace(Nurbs2D.knots.U{k}(iu),Nurbs2D.knots.U{k}(iu+1),np_u);
        v_sample = linspace(Nurbs2D.knots.V{k}(iv),Nurbs2D.knots.V{k}(iv+1),np_v);
        x_sample = zeros(np_u*np_v,1);
        y_sample = zeros(np_u*np_v,1);
        z_sample = zeros(np_u*np_v,1);
        count = 1;
        for j = 1:np_v
            for r = 1:np_u
                dNu = dersbasisfuns(iu,u_sample(r),Nurbs2D.order{k}(1)-1,0,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,0,Nurbs2D.knots.V{k});
                CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                Sw = zeros(4,1);
                for i = 1:4
                    Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{k}(1),Nurbs2D.order{k}(2))*dNv(1,:)';
                end
                S = Sw(1:3,:) / Sw(4);
                x_sample(count) = S(1);
                y_sample(count) = S(2);
                z_sample(count) = S(3);
                elData(count,:,count_el) = [S(1) S(2) S(3)];
                nodeData(count_node,:) = [S(1), S(2), S(3)];
                count = count+1;
                count_node = count_node+1;
            end
        end
        hold on
        scatter3(x_sample,y_sample,z_sample,80,'b','filled');
        hold off
        count_el = count_el+1;
    end
end
axis off
%---------------------------------------------------------------
% Patch elementsectivity
% Nodal Coordinates (nodes)
% elementsectivity (elements)
%---------------------------------------------------------------
TOL = 0.005; %---> Check!
nodes_sem = uniquetol(nodeData,TOL,'ByRows',true);
elements_sem = zeros(tot_el,np_u*np_v);
for i = 1:tot_el
    for j = 1:np_u*np_v
        node_id = find(ismembertol(nodes_sem, elData(j,:,i),TOL,'ByRows',true));
        elements_sem(i,j) = node_id;
    end
end
save('nodes.mat','nodes_sem');
save('elements.mat','elements_sem');
%---------------------------------------
clear
close all
clc
tic
%% -----------------------------------------------------------------------
% ------------------- Input Parameters --------------------------------
% -------------------------------------------------------------------------
% Elements and nodes
load elements
load nodes
elements = elements_sem;
nodes = nodes_sem;
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
%% ------------------------------------------------------------------------
% -------------- Assembly -------------------------------------------------
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
disp(['Assembly: ' num2str(round(toc,1)) ' s'])
%% ------------------------------------------------------------------------
% -------------- Eigenvalue Solution --------------------------------------
% -------------------------------------------------------------------------
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
drawModeShape2;
% % 
% % ind_bc = unique([find(nodes(:,2)<0.0001); find(nodes(:,1)<0.0001); find(abs(1.4-nodes(:,2))<0.0001); find(abs(nodes(:,1)-2.0)<0.0001)]);
% % ind_dof = transpose(setdiff([1:length(nodes)],ind_bc));
% % U_Modes = zeros(5*length(nodes),size(U,2));
% % U_Modes([ind_dof; length(nodes)+ind_dof; 2*length(nodes)+ind_dof; 3*length(nodes)+ind_dof; 4*length(nodes)+ind_dof],:) = UN;
% % % %------------------
% % H = zeros(size(elements,1),size(elements,1));
% % G = zeros(size(elements,1),size(elements,1));
% % C = 0.5.*eye(size(elements,1));
% % modeNum=20;
% % b = zeros(size(elements,1),modeNum);
% % [xgp,wgp,ngp] = gaussQuad2d(8,8);
% % [N, dN] = linear2Dshapefun(xgp(:,1),xgp(:,2));
% % mid = 13;
% % subd_in = [1 11 13 3
% %                    11 21 23 13
% %                    13 23 25 15
% %                    3 13 15 5];
% % subd_out = [1 21 25 5];
% % yfs = 1.4;
% % for i=1:size(elementpoints,1)
% %     disp(i)
% %     node_i=posn0(elementpoints(i,mid),1:2);
% %     node_ip = [node_i(1), 2*yfs-node_i(2)];
% %     midpoint = elementpoints(i,13);
% %     aa = find((indB==midpoint)&(indA==3));
% %     if ~isempty(aa)
% %         b(i,:) = UN(aa,:);
% %     end
% %     for j=1:size(elementpoints,1)
% %         if j==i
% %             for p=1:4
% %                 xn = posn0(elementpoints(j,subd_in(p,:)),1);
% %                 yn = posn0(elementpoints(j,subd_in(p,:)),2);
% %                 for k=1:ngp
% %                     r_vector = [N(k,:)*xn; N(k,:)*yn];
% %                     J_mat = [dN(1,:,k)*xn, dN(1,:,k)*yn; dN(2,:,k)*xn, dN(2,:,k)*yn];
% %                     J = det(J_mat);
% %                     r=norm(r_vector-transpose(node_i));
% %                     rp = norm(r_vector-transpose(node_ip));
% %                     G(i,j)=G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
% %                 end
% %             end
% %         else
% %             xn = posn0(elementpoints(j,subd_out),1);
% %             yn = posn0(elementpoints(j,subd_out),2);
% %             for k=1:ngp
% %                 r_vector = [N(k,:)*xn; N(k,:)*yn];
% %                 J_mat = [dN(1,:,k)*xn, dN(1,:,k)*yn; dN(2,:,k)*xn, dN(2,:,k)*yn];
% %                 J = det(J_mat);
% %                 r=norm(r_vector-transpose(node_i));
% %                 rp = norm(r_vector-transpose(node_ip));
% %                 G(i,j)=G(i,j)+(1/(4*pi*r)-1/(4*pi*rp))*wgp(k)*J;
% %             end
% %         end
% %     end
% % end
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