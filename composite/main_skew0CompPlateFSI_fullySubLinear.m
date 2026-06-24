%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element Method (SEM)
% Meshing       : NURBS-Based Coarse-Quad Meshing
%
% Description   :
%  This code performs hydroelastic vibration analysis of vertical 
%  rectangular composite plates lie on x-y plane and partially in 
%  contact with fluid from one side using coupled SEM-BEM framework 
%  along with a NURBS-based coarse-quad meshing tecnique.
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
FileName = 'skew0_5x5_';
semPatch = [1]; %Enter # SEM Patches
bemPatch = [2]; %Enter # BEM Patches
numPatch = 2;
%-Thickness:
t = 0.01;
psi = 30;
%-Stacking Sequence:
stSeq = [45 -45 45 -45 45];
%-Material ID:
materialID = 'rectComposite';
%-Order of SEM elements:
N = 5;
%-Number of Modes to be extracted:
modeNum = 48;
%-Number of Modes to be plotted:
modeNumPlot = 4;
%-Element Type
ET = 3; % 1: Plate on x-y plane (3 DOF)
        % 2: Shell in 3D (6 DOF)
        % 3: Composite Plate on x-y plane (5 DOF)
        % 4: Composite Shell in 3D (6 DOF)
%-DOF per Sampling Point:
shell_dof = 5;
fluid_dof = 1;
%-Formulation
form = 1;
%----------------
% Pre-Processing
%----------------
tic;
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
toc;
%
tic;
sembem2D = sembem2Dmesh(Nurbs2D,N,semPatch,bemPatch,shell_dof,fluid_dof);
sembem2D.t = t;
sembem2D.stSeq = stSeq;
sembem2D.ET = ET;
%sembem2D.shell_dof = shell_dof;
sembem2D.form = form;
sembem2D.N = N;
sembem2D.non = size(sembem2D.nodes,1);
sembem2D.dof = sembem2D.non*shell_dof;
sembem2D.matProp = get_material_property(materialID);
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
%
ind1 = find(sembem2D.nodes(:,1)<x_min+1E-5|sembem2D.nodes(:,1)>x_max-1E-5);
ind2 = find(sembem2D.nodes(:,2)<y_min+1E-5|sembem2D.nodes(:,2)>y_max-1E-5);
%CCCC
BounNodes = unique([5.*ind1-4; 5.*ind1-3; 5.*ind1-2; 5.*ind1-1; 5.*ind1; ...
                    5.*ind2-4; 5.*ind2-3; 5.*ind2-2; 5.*ind2-1; 5.*ind2]);
%SSSS
% BounNodes = unique([5.*ind1-4; 5.*ind1-3; 5.*ind1-2; 5.*ind2-4; 5.*ind2-3; 5.*ind2-2]);
%
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
toc;
tic;
%-Eigenvalue Solver
sigma = 0.01;
[V,freq] = eigs(K,M,modeNum,sigma);
[freq,loc] = sort((sqrt(diag(freq)-sigma)));
%
toc;
V = V(:,loc);
freqHz = freq/2/pi;
%
U = zeros(size(V,1), size(V,2));
for i=1:size(V,2)
    U(:,i) = V(:,loc(i));
end
% Mass normalization of mode shapes
UN = zeros(size(U,1),size(U,2));
Mnorm = diag(U'*(M)*U).^(1/2);
for i=1:size(U,2)
    UN(:,i)=U(:,i)/Mnorm(i);
end
%
all_nodes = 1:sembem2D.dof;
active = setdiff(all_nodes,BounNodes);
uModes = zeros(sembem2D.dof,modeNum);
uModes(active,1:modeNum) = uModes(active,1:modeNum) + UN(:,1:modeNum);
sembem2D.uModes = uModes;
sembem2D.freq = freq;
sembem2D.freqHz = freqHz;
% % %-----------------
% % % Post-Processing
% % %-----------------
plotModeShapes(sembem2D,modeNumPlot);
%---------------------------------------------
% Displacement Amplitudes at BEM nodes (in Z)
%---------------------------------------------
uModesZ = zeros(size(sembem2D.nodes,1),modeNum);
rModesZ = zeros(size(sembem2D.nodes,1),3);
%
space_xi.a = -1; space_xi.b = 1;
xi_Int = space_xi.a:(space_xi.b-space_xi.a)/(N-1):space_xi.b;
space_xi.N = N;
space_eta.a = -1; space_eta.b = 1;
eta_Int = space_eta.a:(space_eta.b-space_eta.a)/(N-1):space_eta.b;
space_eta.N = N;
%
for i = 1:modeNum
    for el = 1:size(sembem2D.conn,1)
        %
        nconn = sembem2D.conn(el,shell_dof:shell_dof:end)./shell_dof;
        el_conn = sembem2D.conn(el,:);
        %
        deflection_w = sembem2D.uModes(el_conn(3:5:end),i);
        %
        a_w = sembem2D.FT*deflection_w;
        %
        deflection_w_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, a_w);
        %
        if i == 1
            xBEM = sembem2D.FT*sembem2D.nodes(nconn,1);
            xBEM_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, xBEM);
            yBEM = sembem2D.FT*sembem2D.nodes(nconn,2);
            yBEM_Int = Interpol2D(space_xi, space_eta, xi_Int, eta_Int, yBEM);
            rModesZ(nconn,1:2) = [xBEM_Int, yBEM_Int];
        end
        %
        uModesZ(nconn,i) = deflection_w_Int;
    end
end
sembem2D.uModesZ = uModesZ;
%
%
for i = 1:modeNumPlot
    figure
    hold on
    %
    if max(sembem2D.uModesZ(:,i)) > -min(sembem2D.uModesZ(:,i))
        modesign = 1;
        clim1 = max(sembem2D.uModesZ(:,i));
    else
        modesign = -1;
        clim1 = -min(sembem2D.uModesZ(:,i));
    end
    %
    for el = 1:sembem2D.nel
        %
        nconn = sembem2D.conn(el,5:5:end)./5;
        %
        x_el = reshape(rModesZ(nconn,1),sembem2D.N,sembem2D.N);
        y_el = reshape(rModesZ(nconn,2),sembem2D.N,sembem2D.N);
        u_el = reshape(sembem2D.uModesZ(nconn,i),sembem2D.N,sembem2D.N);
        %
        surf(x_el,y_el,modesign.*u_el)
        %
    end
    %
    hold off
    axis equal
    axis off
    box on
    shading flat
    colormap jet
    caxis([-clim1 clim1]);
    view(0,90)
    title(['Mode ' num2str(i)],'FontSize',12,'FontWeight','normal')
end
%
%
% % tic;
% % countBEM = size(sembem2D.nodesBEM,1);
% % %
% % % Bem matrices and vectors
% % H = zeros(countBEM,countBEM);
% % G = zeros(countBEM,countBEM);
% % C = 0.5.*eye(countBEM,countBEM);
% % b = zeros(countBEM,modeNum);
% % %------------------------------------
% % % Gaussian Quadrature
% % [xgp,wgp,ngp] = gaussQuad2d(4,4);
% % %------------------------------------
% % % Tolerance
% % dist_tol = 2;
% % %------------------------------------
% % pBEM = 1;
% % [N, dN] = shapefunc2D(xgp(:,1),xgp(:,2),pBEM);
% % count_col = 1;
% % %
% % for j=1:size(sembem2D.nodesBEM,1)
% %     %
% %     node_i = sembem2D.nodesBEM(j,:);
% %     %
% %     ni = [0,0,1];
% %     ind = find((abs(rModesZ(:,1)-node_i(1))<1E-4) & (abs(rModesZ(:,2)-node_i(2))<1E-4)); %& (abs(rModesZ(:,3)-node_i(3)))<1E-5);
% %     %
% %     if ~isempty(ind)
% %         b(count_col,:) = sembem2D.uModesZ(ind,:);
% %     else
% %         disp('b empty')
% %     end
% %     %
% %     for m=1:size(sembem2D.connBEM,1)
% %         %
% %         xn = sembem2D.nodesBEM(sembem2D.connBEM(m,:),1);
% %         yn = sembem2D.nodesBEM(sembem2D.connBEM(m,:),2);
% %         zn = sembem2D.nodesBEM(sembem2D.connBEM(m,:),3);
% %         %
% %         for g=1:ngp
% %             %
% %             posj = [N(g,:)*xn; N(g,:)*yn; N(g,:)*zn];
% %             a1j = [dN(1,:,g)*xn; dN(1,:,g)*yn; dN(1,:,g)*zn];
% %             a2j = [dN(2,:,g)*xn; dN(2,:,g)*yn; dN(2,:,g)*zn];
% %             J = norm(cross(a1j,a2j));
% %             nj = cross(a1j,a2j)./J;
% %             r_vector = posj-transpose(node_i);
% %             r=norm(posj-transpose(node_i));
% %             dr_dn = (r_vector'*nj)/r;
% %             G(count_col,sembem2D.connBEM(m,:)) = G(count_col,sembem2D.connBEM(m,:)) + ((1/(4*pi*r))*wgp(g)*J).*N(g,:);
% %             H(count_col,sembem2D.connBEM(m,:)) = H(count_col,sembem2D.connBEM(m,:)) + (((-1/(4*pi*r^2))*dr_dn)*wgp(g)*J).*N(g,:);
% %         end
% %     end
% %     count_col = count_col + 1
% % end
% % %
% % phi = (C-H)\(G*b);
% % a = zeros(modeNum,modeNum);
% % for i = 1:modeNum
% %     for el = 1:size(sembem2D.connBEM,1)
% %         %
% %         phi_el = phi(sembem2D.connBEM(el,:),i);
% %         eigvec_el = b(sembem2D.connBEM(el,:),:);
% %         xn = sembem2D.nodesBEM(sembem2D.connBEM(el,:),1);
% %         yn = sembem2D.nodesBEM(sembem2D.connBEM(el,:),2);
% %         zn = sembem2D.nodesBEM(sembem2D.connBEM(el,:),3);
% %         for g=1:ngp
% %             %
% %             posj = [N(g,:)*xn; N(g,:)*yn; N(g,:)*zn];
% %             a1j = [dN(1,:,g)*xn; dN(1,:,g)*yn; dN(1,:,g)*zn];
% %             a2j = [dN(2,:,g)*xn; dN(2,:,g)*yn; dN(2,:,g)*zn];
% %             J = norm(cross(a1j,a2j));
% %             phi_g = N(g,:)*phi_el;
% %             eigvec_g = N(g,:)*eigvec_el;
% %             for j = 1:modeNum
% %                 a(i,j) = a(i,j) + (1000*wgp(g)*J).*(phi_g*(eigvec_g(j)));
% %             end
% %         end
% %     end
% % end
% % A = eye(modeNum,modeNum);
% % c = zeros(modeNum,modeNum);
% % for i = 1:modeNum
% %     c(i,i) = sembem2D.freq(i)^2;
% % end
% % [wV, wfreq] = eig(c,(a+A));
% % wfreq = diag(wfreq);
% % [wfreq2,ind] = sort((sqrt(real(wfreq))./(2*pi)));
% % wetV = wV(:,ind);
% % toc;
% % tEnd = cputime-tStart;
% % % ------------------------------------------------------------------------
% % % ----- PLOT WET MODES ---------------------------------------------------
% % % ------------------------------------------------------------------------
% % numWmode = 8;
% % for k = 1:numWmode
% %     wetModeDisp = 0.*sembem2D.uModes(:,k);
% %     for j = 1:modeNum
% %         wetModeDisp = wetModeDisp - (wetV(j,k)).*sembem2D.uModes(:,j);
% %     end
% %     if max(wetModeDisp(3:5:end)) > -min(wetModeDisp(3:5:end))
% %         modesign = 1;
% %         clim1 = max(wetModeDisp(3:5:end));
% %     else
% %         modesign = -1;
% %         clim1 = -min(wetModeDisp(3:5:end));
% %     end
% %     % clim1 = max(abs(min(wetModeDisp(indW))),abs(max(wetModeDisp(indW))));
% %     figure;
% %     hold on
% %     for el = 1:sembem2D.nel
% %         %
% %         nconn = sembem2D.conn(el,5:5:end)./5;
% %         %
% %         x_el = reshape(sembem2D.nodes(nconn,1),sembem2D.N,sembem2D.N);
% %         y_el = reshape(sembem2D.nodes(nconn,2),sembem2D.N,sembem2D.N);
% %         u_el = reshape(wetModeDisp(sembem2D.conn(el,3:5:end)),sembem2D.N,sembem2D.N);
% %         %
% %         surf(x_el,y_el,modesign.*u_el)
% %         %
% %     end
% %     hold off
% %     axis equal
% %     axis off
% %     box on
% %     shading interp
% %     colormap jet
% %     caxis([-clim1 clim1])
% %     view(0,90)
% % end