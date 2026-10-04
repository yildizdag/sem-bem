function sembem2D = sembem2DmeshAdaptive(Nurbs2D,delta0,a_s,...
                                         semPatch,bemPatch,...
                                         shell_dof,fluid_dof)

%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Chebyshev Spectral Element Method (SEM)
% Meshing       : NURBS-Based Adaptive Coarse-Quad Meshing
%
% Description   :
%  This function generates a NURBS-based SEM-BEM mesh. The number of
%  spectral sampling points is determined adaptively using a geometric
%  convergence criterion.
%
%  For each NURBS element, the procedure starts from N = 3 sampling points
%  in each parametric direction. The physical distances between consecutive
%  sampling points along the four element edges are calculated.
%
%  For each edge:
%
%       a_e       = sum(d_i)
%       Delta_max = max(d_i)
%       delta     = Delta_max / sqrt(a_s*a_e)
%
%  N is increased until the convergence criterion is satisfied for all
%  four edges of the element. After determining the required N for every
%  element, the maximum value is used as the common SEM discretization.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================


%% General information ----------------------------------------------------

numSEMpatch = size(semPatch,2);
numBEMpatch = size(bemPatch,2);

nel = 0;

for k = 1:numSEMpatch
    nel = nel + Nurbs2D.nel{k};
end

sembem2D.nel       = nel;
sembem2D.shell_dof = shell_dof;
sembem2D.fluid_dof = fluid_dof;


%% Adaptive determination of sampling points ------------------------------

% Minimum number of sampling points per direction
N_min = 3;

% Safety limit to avoid an infinite loop in case of an inappropriate
% convergence threshold or problematic geometry
N_max = 50;

% Storage for element-wise adaptive results
N_element     = zeros(nel,1);
delta_element = zeros(nel,1);

% Storage of directional/edge information at convergence
delta_edges = zeros(nel,4);
a_e_edges   = zeros(nel,4);

count_el = 1;


for k = 1:numSEMpatch

    for el = 1:Nurbs2D.nel{k}

        %% Parametric domain of current NURBS element ---------------------

        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);

        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);

        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);


        %% Control points of current NURBS element ------------------------

        CP = Nurbs2D.cPoints{k}(:,...
             iu-Nurbs2D.order{k}(1)+1:iu,...
             iv-Nurbs2D.order{k}(2)+1:iv);


        %% Start adaptive iteration --------------------------------------

        N_test = N_min;

        while true

            %% Lobatto points --------------------------------------------

            xi_test  = lobat(N_test);
            eta_test = lobat(N_test);


            %% Mapping to NURBS parametric domain ------------------------

            u_sample_test = ...
                0.5.*(1-xi_test).*u1 + ...
                0.5.*(1+xi_test).*u2;

            v_sample_test = ...
                0.5.*(1-eta_test).*v1 + ...
                0.5.*(1+eta_test).*v2;


            %% Physical coordinates of N x N sampling points ------------

            samplePoints = zeros(N_test,N_test,3);

            for i = 1:N_test

                for j = 1:N_test

                    dNu = dersbasisfuns(...
                          iu,...
                          u_sample_test(i),...
                          Nurbs2D.order{k}(1)-1,...
                          2,...
                          Nurbs2D.knots.U{k});

                    dNv = dersbasisfuns(...
                          iv,...
                          v_sample_test(j),...
                          Nurbs2D.order{k}(2)-1,...
                          2,...
                          Nurbs2D.knots.V{k});

                    [~,dS] = derRat2DBasisFuns(...
                              dNu,...
                              dNv,...
                              Nurbs2D.order{k}(1),...
                              Nurbs2D.order{k}(2),...
                              CP,...
                              2,...
                              2);

                    samplePoints(i,j,:) = dS(:,1,1);

                end

            end


            %% Extract four physical element edges -----------------------

            % Edge 1: v = v1
            edge1 = squeeze(samplePoints(:,1,:));

            % Edge 2: u = u2
            edge2 = squeeze(samplePoints(end,:,:));

            % Edge 3: v = v2
            edge3 = squeeze(samplePoints(:,end,:));

            % Edge 4: u = u1
            edge4 = squeeze(samplePoints(1,:,:));

            edges = {edge1,edge2,edge3,edge4};


            %% Calculate convergence parameter for four edges ------------

            a_e_now = zeros(4,1);
            delta_max_now = zeros(4,1);
            delta_now = zeros(4,1);

            for edge = 1:4

                % Physical distances between consecutive sampling points
                d = vecnorm(diff(edges{edge},1,1),2,2);

                % Approximate physical length of the NURBS element edge
                a_e_now(edge) = sum(d);

                % Maximum distance between consecutive sampling points
                delta_max_now(edge) = max(d);

                % Dimensionless convergence parameter
                delta_now(edge) = ...
                    delta_max_now(edge) / ...
                    sqrt(a_s*a_e_now(edge));

            end


            %% Governing delta for current element -----------------------

            delta_el = max(delta_now);


            %% Check convergence -----------------------------------------

            if delta_el <= delta0

                break

            end


            %% Increase sampling resolution ------------------------------

            N_test = N_test + 1;


            %% Safety check ----------------------------------------------

            if N_test > N_max

                error(['Adaptive SEM sampling did not converge. ',...
                       'N exceeded N_max = %d for patch %d, element %d.'],...
                       N_max,k,el);

            end

        end


        %% Store converged element information ---------------------------

        N_element(count_el) = N_test;

        delta_element(count_el) = delta_el;

        delta_edges(count_el,:) = delta_now.';

        a_e_edges(count_el,:) = a_e_now.';

        count_el = count_el + 1;

    end

end


%% Determine common SEM sampling number ----------------------------------

N = max(N_element);

sembem2D.N = N;

% Polynomial degree corresponding to N Lobatto sampling points
sembem2D.p = N-1;

% Store adaptive information for post-processing/reporting
sembem2D.N_element     = N_element;
sembem2D.delta_element = delta_element;
sembem2D.delta_edges   = delta_edges;
sembem2D.a_e_edges     = a_e_edges;

sembem2D.delta0 = delta0;
sembem2D.a_s    = a_s;


%% Display adaptive results -----------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ADAPTIVE SEM SAMPLING\n');
fprintf('============================================================\n');
fprintf('delta0                  = %.6f\n',delta0);
fprintf('Structural length a_s   = %.6f\n',a_s);
fprintf('Minimum sampling N_min  = %d\n',N_min);
fprintf('------------------------------------------------------------\n');
fprintf('Element        N_e          p_e        delta_e\n');
fprintf('------------------------------------------------------------\n');

for el = 1:nel

    fprintf('%5d          %3d          %3d       %.6f\n',...
            el,...
            N_element(el),...
            N_element(el)-1,...
            delta_element(el));

end

fprintf('------------------------------------------------------------\n');
fprintf('Global sampling number N = %d\n',N);
fprintf('Global polynomial degree p = %d\n',N-1);
fprintf('============================================================\n\n');


%% SEM mesh generation ----------------------------------------------------

ntot = N*N*nel;

nodeData = zeros(ntot,3);

JacMatData = zeros(3,2,N*N,nel);
InvJacMatData = zeros(2,2,N*N,nel);
JacobianData = zeros(1,N*N,nel);

curvData = zeros(ntot,2);

T1Data = zeros(ntot,3);
T2Data = zeros(ntot,3);
NData  = zeros(ntot,3);

RData = zeros(3,3,ntot);

count_el   = 1;
count_node = 1;

xi  = lobat(N);
eta = lobat(N);

epsilon = 1E-6;


for k = 1:numSEMpatch

    for el = 1:Nurbs2D.nel{k}

        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);

        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);

        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);

        u_sample = ...
            0.5.*(1-xi).*u1 + ...
            0.5.*(1+xi).*u2;

        v_sample = ...
            0.5.*(1-eta).*v1 + ...
            0.5.*(1+eta).*v2;

        count = 1;

        CP = Nurbs2D.cPoints{k}(:,...
             iu-Nurbs2D.order{k}(1)+1:iu,...
             iv-Nurbs2D.order{k}(2)+1:iv);

        du = ...
            (Nurbs2D.knots.U{k}(iu+1) - ...
             Nurbs2D.knots.U{k}(iu))/2;

        dv = ...
            (Nurbs2D.knots.V{k}(iv+1) - ...
             Nurbs2D.knots.V{k}(iv))/2;


        for i = 1:N

            for j = 1:N

                dNu = dersbasisfuns(...
                      iu,...
                      u_sample(i),...
                      Nurbs2D.order{k}(1)-1,...
                      2,...
                      Nurbs2D.knots.U{k});

                dNv = dersbasisfuns(...
                      iv,...
                      v_sample(j),...
                      Nurbs2D.order{k}(2)-1,...
                      2,...
                      Nurbs2D.knots.V{k});

                [~,dS] = derRat2DBasisFuns(...
                          dNu,...
                          dNv,...
                          Nurbs2D.order{k}(1),...
                          Nurbs2D.order{k}(2),...
                          CP,...
                          2,...
                          2);

                nodeData(count_node,:) = ...
                    epsilon.*(dS(:,1,1)'./epsilon);


                %% Surface basis vectors ---------------------------------

                A1 = dS(:,2,1);
                A2 = dS(:,1,2);

                t1 = A1 ./ norm(A1);

                A3 = cross(A1,A2) / ...
                     norm(cross(A1,A2));

                t2 = cross(A3,t1);
                t2 = t2 ./ norm(t2);


                %% Curvature calculation ---------------------------------

                F1 = [A1 A2]'*[A1 A2];

                Ac = [A1, A2]/F1;

                F2 = ...
                    [dot(dS(:,3,1),A3), ...
                     dot(dS(:,2,2),A3); ...
                     dot(dS(:,2,2),A3), ...
                     dot(dS(:,1,3),A3)];

                F = F1\F2;

                kappa = eig(F);


                %% Store geometrical quantities ---------------------------

                JacMatData(:,:,count,count_el) = [A1, A2];

                JacobianData(1,count,count_el) = ...
                    norm(cross(A1,A2))*du*dv;

                InvJacMatData(:,:,count,count_el) = ...
                    [dot(t1,Ac(:,1))/du, ...
                     dot(t1,Ac(:,2))/dv; ...
                     dot(t2,Ac(:,1))/du, ...
                     dot(t2,Ac(:,2))/dv];

                curvData(count_node,:) = ...
                    [abs(kappa(1)),abs(kappa(2))];

                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';

                RData(:,:,count_node) = ...
                    [t1,t2,A3];


                count = count + 1;

                count_node = count_node + 1;

            end

        end

        count_el = count_el + 1;

    end

end


%% Merge coincident SEM nodes ---------------------------------------------

TOL = 1e-5;

[nodes_sem,IA,IC] = ...
    uniquetol(nodeData,TOL,'ByRows',true);

Kappa = curvData(IA,:);

elemNode = reshape(IC,N*N,nel).';

conn_sem = zeros(nel,shell_dof*N*N);

for d = 1:shell_dof

    conn_sem(:,d:shell_dof:end) = ...
        shell_dof*elemNode-(shell_dof-d);

end


%% Store SEM mesh data ----------------------------------------------------

sembem2D.nodes = nodes_sem;

sembem2D.conn = conn_sem;

sembem2D.Jmat = JacMatData;

sembem2D.J = JacobianData;

sembem2D.InvJmat = InvJacMatData;

sembem2D.Kappa = Kappa;

sembem2D.t1 = T1Data(IA,:);

sembem2D.t2 = T2Data(IA,:);

sembem2D.n = NData(IA,:);

sembem2D.R = RData(:,:,IA);


%% Chebyshev operators ----------------------------------------------------

space.a = -1;
space.b = 1;
space.N = N;

[FT_xi,BT_xi] = cheb(space);

D_xi = derivative(space);

V_xi = InnerProduct(space);

Q1_xi = BT_xi*D_xi*FT_xi;

Q2_xi = BT_xi*D_xi^2*FT_xi;


[FT_eta,BT_eta] = cheb(space);

D_eta = derivative(space);

V_eta = InnerProduct(space);

Q1_eta = BT_eta*D_eta*FT_eta;

Q2_eta = BT_eta*D_eta^2*FT_eta;


sembem2D.VD = kron(V_xi,V_eta);

sembem2D.Q1xi = ...
    kron(Q1_xi,eye(N));

sembem2D.Q1eta = ...
    kron(eye(N),Q1_eta);

sembem2D.Q2xi = ...
    kron(Q2_xi,eye(N));

sembem2D.Q2eta = ...
    kron(eye(N),Q2_eta);

sembem2D.Qxieta = ...
    sembem2D.Q2xi*sembem2D.Q2eta;

sembem2D.FT = ...
    Fxy_mapping(N,N,FT_xi,FT_eta);


%% BEM mesh generation ----------------------------------------------------
%
% IMPORTANT:
% Keep the original BEM section from sembem2Dmesh.m below this point.
%
% The adaptive procedure above determines the common SEM sampling number N.
% No modification to the existing BEM formulation is required at this
% stage.
%

if numBEMpatch > 0

    N_BEM = 5;

    % -------------------------------------------------------------
    % PASTE THE ORIGINAL BEM SECTION HERE WITHOUT MODIFICATION
    % -------------------------------------------------------------

end


end