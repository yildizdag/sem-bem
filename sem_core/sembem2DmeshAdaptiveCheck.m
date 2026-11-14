function sembem2D = sembem2DmeshAdaptiveCheck(Nurbs2D,delta0,a_s,...
                                               semPatch,bemPatch,...
                                               shell_dof,fluid_dof)

%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
% Method        : Hybrid SEM-BEM
% Meshing       : NURBS-Based Adaptive Coarse-Quad Meshing
%
% Description   :
%  This function generates a NURBS-based SEM-BEM mesh. The number of
%  SEM-BEM sampling points is determined adaptively using a geometric
%  convergence criterion.
%
%  In addition, two geometrical checks are performed:
%
%    1) Jacobian positivity / degeneracy check
%    2) Surface-normal orientation consistency check
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


%% Adaptive determination of SEM-BEM sampling points ----------------------

N_min = 3;
N_max = 50;

N_element     = zeros(nel,1);
delta_element = zeros(nel,1);

delta_edges = zeros(nel,4);
a_e_edges   = zeros(nel,4);

count_el = 1;

for k = 1:numSEMpatch

    for el = 1:Nurbs2D.nel{k}

        %% Parametric domain ---------------------------------------------

        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);

        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);

        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);


        %% Control points ------------------------------------------------

        CP = Nurbs2D.cPoints{k}(:,...
             iu-Nurbs2D.order{k}(1)+1:iu,...
             iv-Nurbs2D.order{k}(2)+1:iv);


        %% Adaptive iteration -------------------------------------------

        N_test = N_min;

        while true

            %% Lobatto points -------------------------------------------

            xi_test  = lobat(N_test);
            eta_test = lobat(N_test);


            %% Mapping to NURBS parametric domain -----------------------

            u_sample_test = ...
                0.5.*(1-xi_test).*u1 + ...
                0.5.*(1+xi_test).*u2;

            v_sample_test = ...
                0.5.*(1-eta_test).*v1 + ...
                0.5.*(1+eta_test).*v2;


            %% Physical sampling points ---------------------------------

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


            %% Four physical element edges ------------------------------

            edge1 = squeeze(samplePoints(:,1,:));
            edge2 = squeeze(samplePoints(end,:,:));
            edge3 = squeeze(samplePoints(:,end,:));
            edge4 = squeeze(samplePoints(1,:,:));

            edges = {edge1,edge2,edge3,edge4};


            %% Convergence parameter ------------------------------------

            a_e_now       = zeros(4,1);
            delta_max_now = zeros(4,1);
            delta_now     = zeros(4,1);

            for edge = 1:4

                d = vecnorm(diff(edges{edge},1,1),2,2);

                a_e_now(edge) = sum(d);

                delta_max_now(edge) = max(d);

                delta_now(edge) = ...
                    delta_max_now(edge) / ...
                    sqrt(a_s*a_e_now(edge));

            end


            %% Governing value ------------------------------------------

            delta_el = max(delta_now);


            %% Convergence check ----------------------------------------

            if delta_el <= delta0
                break
            end

            N_test = N_test + 1;


            %% Safety check ---------------------------------------------

            if N_test > N_max

                error(['Adaptive SEM-BEM sampling did not converge. ',...
                       'N exceeded N_max = %d for patch %d, element %d.'],...
                       N_max,k,el);

            end

        end


        %% Store adaptive information -----------------------------------

        N_element(count_el)     = N_test;
        delta_element(count_el) = delta_el;

        delta_edges(count_el,:) = delta_now.';
        a_e_edges(count_el,:)   = a_e_now.';

        count_el = count_el + 1;

    end

end


%% Common SEM-BEM sampling number -----------------------------------------

N = max(N_element);

sembem2D.N = N;
sembem2D.p = N-1;

sembem2D.N_element     = N_element;
sembem2D.delta_element = delta_element;
sembem2D.delta_edges   = delta_edges;
sembem2D.a_e_edges     = a_e_edges;

sembem2D.delta0 = delta0;
sembem2D.a_s    = a_s;


%% Display adaptive results -----------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ADAPTIVE SEM-BEM SAMPLING\n');
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


%% SEM-BEM mesh generation ------------------------------------------------

ntot = N*N*nel;

nodeData = zeros(ntot,3);

JacMatData    = zeros(3,2,N*N,nel);
InvJacMatData = zeros(2,2,N*N,nel);
JacobianData  = zeros(1,N*N,nel);

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

        %% Parametric element domain ------------------------------------

        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);

        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);

        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);


        %% SEM-BEM -> NURBS parametric mapping --------------------------

        u_sample = ...
            0.5.*(1-xi).*u1 + ...
            0.5.*(1+xi).*u2;

        v_sample = ...
            0.5.*(1-eta).*v1 + ...
            0.5.*(1+eta).*v2;


        %% Control points ------------------------------------------------

        CP = Nurbs2D.cPoints{k}(:,...
             iu-Nurbs2D.order{k}(1)+1:iu,...
             iv-Nurbs2D.order{k}(2)+1:iv);


        %% Parametric transformation ------------------------------------

        du = (u2-u1)/2;
        dv = (v2-v1)/2;

        count = 1;


        %% SEM-BEM sampling points --------------------------------------

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


                %% Physical point ---------------------------------------

                nodeData(count_node,:) = ...
                    epsilon.*(dS(:,1,1)'./epsilon);


                %% Surface basis vectors --------------------------------

                A1 = dS(:,2,1);
                A2 = dS(:,1,2);

                crossA    = cross(A1,A2);
                crossNorm = norm(crossA);

                if crossNorm > eps

                    t1 = A1 ./ norm(A1);

                    A3 = crossA ./ crossNorm;

                    t2 = cross(A3,t1);
                    t2 = t2 ./ norm(t2);

                else

                    t1 = [NaN;NaN;NaN];
                    t2 = [NaN;NaN;NaN];
                    A3 = [NaN;NaN;NaN];

                end


                %% Curvature calculation --------------------------------

                if crossNorm > eps

                    F1 = [A1 A2]'*[A1 A2];

                    Ac = [A1,A2]/F1;

                    F2 = ...
                        [dot(dS(:,3,1),A3), ...
                         dot(dS(:,2,2),A3); ...
                         dot(dS(:,2,2),A3), ...
                         dot(dS(:,1,3),A3)];

                    F = F1\F2;

                    kappa = eig(F);

                else

                    Ac    = zeros(3,2);
                    kappa = [NaN;NaN];

                end


                %% Store geometrical quantities -------------------------

                JacMatData(:,:,count,count_el) = [A1,A2];

                JacobianData(1,count,count_el) = ...
                    crossNorm*du*dv;

                if crossNorm > eps

                    InvJacMatData(:,:,count,count_el) = ...
                        [dot(t1,Ac(:,1))/du, ...
                         dot(t1,Ac(:,2))/dv; ...
                         dot(t2,Ac(:,1))/du, ...
                         dot(t2,Ac(:,2))/dv];

                else

                    InvJacMatData(:,:,count,count_el) = NaN(2,2);

                end

                curvData(count_node,:) = ...
                    [abs(kappa(1)),abs(kappa(2))];

                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';

                RData(:,:,count_node) = [t1,t2,A3];

                count      = count + 1;
                count_node = count_node + 1;

            end

        end

        count_el = count_el + 1;

    end

end


%% Jacobian check ---------------------------------------------------------

jacobianTolerance = 1E-12;

JacobianMin   = zeros(nel,1);
JacobianValid = true(nel,1);

for el = 1:nel

    J_el = squeeze(JacobianData(1,:,el));

    JacobianMin(el) = min(J_el);

    if ~isfinite(JacobianMin(el)) || ...
            JacobianMin(el) <= jacobianTolerance

        JacobianValid(el) = false;

    end

end

invalidJacobianElements = find(~JacobianValid);

meshJacobianValid = isempty(invalidJacobianElements);


%% Merge coincident SEM-BEM nodes ----------------------------------------

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


%% Element-wise normal data -----------------------------------------------

NormalElement = zeros(N*N,3,nel);

for el = 1:nel

    firstID = (el-1)*N*N + 1;
    lastID  = el*N*N;

    NormalElement(:,:,el) = ...
        NData(firstID:lastID,:);

end


%% Surface orientation connectivity ---------------------------------------

%--------------------------------------------------------------------------
% Jacobian magnitude:
%
%       J = ||A1 x A2|| du dv
%
% cannot distinguish two opposite surface orientations.
%
% Therefore, neighboring elements are connected through their common
% SEM-BEM edge nodes. The average dot product of the surface normals on
% each common edge defines the relative orientation:
%
%       +1 : same orientation
%       -1 : opposite orientation
%
% These relative orientations are then propagated through the mesh.
%--------------------------------------------------------------------------

orientationAdjacency = zeros(nel,nel);

MinNormalDot = ones(nel,1);

for elA = 1:nel-1

    nodesA = elemNode(elA,:);

    for elB = elA+1:nel

        nodesB = elemNode(elB,:);

        [commonNodes,idxA,idxB] = ...
            intersect(nodesA,nodesB); %#ok<ASGLU>


        % Elements sharing only a corner are not edge neighbors.
        % A complete SEM-BEM edge contains N nodes.

        if length(commonNodes) < N
            continue
        end


        %% Compare normals on common edge --------------------------------

        edgeDots = zeros(length(commonNodes),1);

        validComparison = true;

        for q = 1:length(commonNodes)

            nA = NormalElement(idxA(q),:,elA);
            nB = NormalElement(idxB(q),:,elB);

            if any(~isfinite(nA)) || any(~isfinite(nB))

                validComparison = false;
                break

            end

            edgeDots(q) = dot(nA,nB);

        end


        if ~validComparison
            continue
        end


        meanDot = mean(edgeDots);

        MinNormalDot(elA) = ...
            min(MinNormalDot(elA),min(edgeDots));

        MinNormalDot(elB) = ...
            min(MinNormalDot(elB),min(edgeDots));


        %% Relative orientation -----------------------------------------

        if meanDot >= 0

            relativeOrientation = 1;

        else

            relativeOrientation = -1;

        end

        orientationAdjacency(elA,elB) = relativeOrientation;
        orientationAdjacency(elB,elA) = relativeOrientation;

    end

end


%% Propagate orientation through connected mesh ---------------------------

% orientationSign:
%
%    +1 -> same orientation as reference element
%    -1 -> reversed relative to reference element
%     0 -> not yet visited

orientationSign = zeros(nel,1);

componentID = zeros(nel,1);

componentCounter = 0;


for startElement = 1:nel

    if orientationSign(startElement) ~= 0
        continue
    end

    componentCounter = componentCounter + 1;

    % First element of each connected component is the reference
    orientationSign(startElement) = 1;
    componentID(startElement) = componentCounter;

    queue = startElement;

    while ~isempty(queue)

        elA = queue(1);
        queue(1) = [];

        neighbors = ...
            find(orientationAdjacency(elA,:) ~= 0);

        for elB = neighbors

            relativeOrientation = ...
                orientationAdjacency(elA,elB);

            expectedSign = ...
                orientationSign(elA)*relativeOrientation;

            if orientationSign(elB) == 0

                orientationSign(elB) = expectedSign;

                componentID(elB) = componentCounter;

                queue(end+1) = elB; %#ok<AGROW>

            end

        end

    end

end


%% Identify reversed elements ---------------------------------------------

reversedElements = find(orientationSign == -1);

OrientationValid = orientationSign == 1;

meshOrientationValid = isempty(reversedElements);


%% Overall mapping status -------------------------------------------------

meshMappingValid = ...
    meshJacobianValid && meshOrientationValid;


%% Display Jacobian check -------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf(' SEM-BEM JACOBIAN CHECK\n');
fprintf('============================================================\n');

fprintf('Element          Minimum Jacobian          Status\n');
fprintf('------------------------------------------------------------\n');

for el = 1:nel

    if JacobianValid(el)
        status = 'VALID';
    else
        status = 'INVALID';
    end

    fprintf('%5d           %14.6e          %s\n',...
            el,...
            JacobianMin(el),...
            status);

end

fprintf('------------------------------------------------------------\n');

if meshJacobianValid

    fprintf('All element Jacobians are positive.\n');

else

    fprintf('Invalid Jacobian element(s): ');

    fprintf('%d ',invalidJacobianElements);

    fprintf('\n');

end

fprintf('============================================================\n\n');


%% Display orientation check ----------------------------------------------

fprintf('============================================================\n');
fprintf(' SEM-BEM GEOMETRIC MAPPING CHECK\n');
fprintf('============================================================\n');

fprintf('Minimum Jacobian          = %.6e\n',...
        min(JacobianMin));


if meshJacobianValid

    fprintf('Jacobian positivity       : PASSED\n');

else

    fprintf('Jacobian positivity       : FAILED\n');

end


if meshOrientationValid

    fprintf('Surface orientation       : PASSED\n');

else

    fprintf('Surface orientation       : FAILED\n');

end

fprintf('------------------------------------------------------------\n');


if meshMappingValid

    fprintf('GEOMETRIC MAPPING CHECKED\n');

else

    if ~meshJacobianValid

        fprintf('Invalid Jacobian element(s): ');

        fprintf('%d ',invalidJacobianElements);

        fprintf('\n');

    end


    if ~meshOrientationValid

        fprintf('Reversed element(s): ');

        fprintf('%d ',reversedElements);

        fprintf('\n');

    end

end

fprintf('============================================================\n\n');


%% Store SEM-BEM mesh data ------------------------------------------------

sembem2D.nodes = nodes_sem;
sembem2D.conn  = conn_sem;

sembem2D.elemNode = elemNode;


%% Jacobian data -----------------------------------------------------------

sembem2D.Jmat    = JacMatData;
sembem2D.J       = JacobianData;
sembem2D.InvJmat = InvJacMatData;

sembem2D.J_min = JacobianMin;

sembem2D.JacobianValid = ...
    JacobianValid;

sembem2D.invalidJacobianElements = ...
    invalidJacobianElements;

sembem2D.meshJacobianValid = ...
    meshJacobianValid;

sembem2D.jacobianTolerance = ...
    jacobianTolerance;


%% Orientation data -------------------------------------------------------

sembem2D.orientationSign = ...
    orientationSign;

sembem2D.orientationValid = ...
    OrientationValid;

sembem2D.reversedElements = ...
    reversedElements;

sembem2D.invalidOrientationElements = ...
    reversedElements;

sembem2D.meshOrientationValid = ...
    meshOrientationValid;

sembem2D.minNormalDot = ...
    MinNormalDot;

sembem2D.orientationAdjacency = ...
    orientationAdjacency;

sembem2D.componentID = ...
    componentID;


%% Overall mapping status -------------------------------------------------

sembem2D.meshMappingValid = ...
    meshMappingValid;


%% Geometric data ---------------------------------------------------------

sembem2D.Kappa = Kappa;

sembem2D.t1 = T1Data(IA,:);
sembem2D.t2 = T2Data(IA,:);
sembem2D.n  = NData(IA,:);

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

% Keep the original BEM implementation here when BEM patches are used.
% The common SEM-BEM sampling number is N.

if numBEMpatch > 0

    N_BEM = N; %#ok<NASGU>

    % Original BEM implementation

end


end