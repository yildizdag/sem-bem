function semOpt2D = semOpt2Dmesh_v1(Nurbs2D_plate,Nurbs2D_stiff,N,shell_dof)
%SEMOPT2DMESH_V1 SEM mesh generation for an optimized plate-stiffener model.
%
% The plate patches are processed first and the stiffener patches second.
% The geometric treatment follows sem2Dmesh_v1.m:
%   - signed principal curvatures,
%   - full signed curvature tensor in the local orthonormal basis,
%   - element-wise local rotation matrices,
%   - metric-degeneracy check.
%

%% Count plate and stiffener elements separately
nel_plate = 0;
for k = 1:Nurbs2D_plate.numpatch
    nel_plate = nel_plate + Nurbs2D_plate.nel{k};
end

nel_stiff = 0;
for k = 1:Nurbs2D_stiff.numpatch
    nel_stiff = nel_stiff + Nurbs2D_stiff.nel{k};
end

nel = nel_plate + nel_stiff;

semOpt2D.nel = nel;
semOpt2D.nel_plate = nel_plate;
semOpt2D.nel_stiff = nel_stiff;
semOpt2D.plateElements = 1:nel_plate;
semOpt2D.stiffElements = nel_plate + (1:nel_stiff);

semOpt2D.N = N;
semOpt2D.shell_dof = shell_dof;

%% Preallocation
ntot = N*N*nel;

nodeData = zeros(ntot,3);
JacMatData = zeros(3,2,N*N,nel);
InvJacMatData = zeros(2,2,N*N,nel);
JacobianData = zeros(1,N*N,nel);

curvData = zeros(2,N*N,nel);
curvTensorData = zeros(2,2,N*N,nel);

T1Data = zeros(ntot,3);
T2Data = zeros(ntot,3);
NData  = zeros(ntot,3);
RData  = zeros(3,3,N*N,nel);

count_el = 1;
count_node = 1;

xi = lobat(N);
eta = lobat(N);

epsilon = 1E-6;

%% ------------------------------------------------------------------------
% Plate patches
% -------------------------------------------------------------------------
for k = 1:Nurbs2D_plate.numpatch

    for el = 1:Nurbs2D_plate.nel{k}

        iu = Nurbs2D_plate.INC{k}(Nurbs2D_plate.IEN{k}(1,el),1);
        iv = Nurbs2D_plate.INC{k}(Nurbs2D_plate.IEN{k}(1,el),2);

        u1 = Nurbs2D_plate.knots.U{k}(iu);
        u2 = Nurbs2D_plate.knots.U{k}(iu+1);
        v1 = Nurbs2D_plate.knots.V{k}(iv);
        v2 = Nurbs2D_plate.knots.V{k}(iv+1);

        u_sample = 0.5.*(1-xi).*u1 + 0.5.*(1+xi).*u2;
        v_sample = 0.5.*(1-eta).*v1 + 0.5.*(1+eta).*v2;

        count = 1;

        CP = Nurbs2D_plate.cPoints{k}( ...
            :, ...
            iu-Nurbs2D_plate.order{k}(1)+1:iu, ...
            iv-Nurbs2D_plate.order{k}(2)+1:iv);

        du = (u2-u1)/2;
        dv = (v2-v1)/2;

        for i = 1:N
            for j = 1:N

                dNu = dersbasisfuns( ...
                    iu,u_sample(i),Nurbs2D_plate.order{k}(1)-1,2, ...
                    Nurbs2D_plate.knots.U{k});

                dNv = dersbasisfuns( ...
                    iv,v_sample(j),Nurbs2D_plate.order{k}(2)-1,2, ...
                    Nurbs2D_plate.knots.V{k});

                [~,dS] = derRat2DBasisFuns( ...
                    dNu,dNv, ...
                    Nurbs2D_plate.order{k}(1), ...
                    Nurbs2D_plate.order{k}(2), ...
                    CP,2,2);

                nodeData(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);

                A1 = dS(:,2,1);
                A2 = dS(:,1,2);

                t1 = A1./norm(A1);

                A3 = cross(A1,A2);
                A3 = A3./norm(A3);

                t2 = cross(A3,t1);
                t2 = t2./norm(t2);

                % First fundamental form
                A = [A1,A2];
                g = A.'*A;

                if rcond(g) < 1E-12
                    error('semOpt2Dmesh_v1:DegeneratePlateMetric', ...
                        ['Degenerate plate surface metric at global element ' ...
                         '%d, sampling point %d.'],count_el,count);
                end

                % Contravariant basis
                Ac = A/g;

                % Signed second fundamental form in parametric coordinates
                A11 = dS(:,3,1);
                A12 = dS(:,2,2);
                A22 = dS(:,1,3);

                bParam = [dot(A11,A3), dot(A12,A3); ...
                          dot(A12,A3), dot(A22,A3)];

                % Curvature tensor in the local orthonormal basis {t1,t2}
                T = [t1,t2];
                C = g\(A.'*T);

                bLocal = C.'*bParam*C;
                bLocal = 0.5.*(bLocal+bLocal.');

                % Signed principal curvatures
                kappa = eig(bLocal);
                kappa = sort(real(kappa));

                JacMatData(:,:,count,count_el) = A;

                JacobianData(1,count,count_el) = ...
                    norm(cross(A1,A2))*du*dv;

                InvJacMatData(:,:,count,count_el) = ...
                    [dot(t1,Ac(:,1))/du, dot(t1,Ac(:,2))/dv; ...
                     dot(t2,Ac(:,1))/du, dot(t2,Ac(:,2))/dv];

                curvData(:,count,count_el) = kappa;
                curvTensorData(:,:,count,count_el) = bLocal;

                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';

                RData(:,:,count,count_el) = [t1,t2,A3];

                count = count + 1;
                count_node = count_node + 1;
            end
        end

        count_el = count_el + 1;
    end
end

%% ------------------------------------------------------------------------
% Stiffener patches
% -------------------------------------------------------------------------
for k = 1:Nurbs2D_stiff.numpatch

    for el = 1:Nurbs2D_stiff.nel{k}

        iu = Nurbs2D_stiff.INC{k}(Nurbs2D_stiff.IEN{k}(1,el),1);
        iv = Nurbs2D_stiff.INC{k}(Nurbs2D_stiff.IEN{k}(1,el),2);

        u1 = Nurbs2D_stiff.knots.U{k}(iu);
        u2 = Nurbs2D_stiff.knots.U{k}(iu+1);
        v1 = Nurbs2D_stiff.knots.V{k}(iv);
        v2 = Nurbs2D_stiff.knots.V{k}(iv+1);

        u_sample = 0.5.*(1-xi).*u1 + 0.5.*(1+xi).*u2;
        v_sample = 0.5.*(1-eta).*v1 + 0.5.*(1+eta).*v2;

        count = 1;

        CP = Nurbs2D_stiff.cPoints{k}( ...
            :, ...
            iu-Nurbs2D_stiff.order{k}(1)+1:iu, ...
            iv-Nurbs2D_stiff.order{k}(2)+1:iv);

        du = (u2-u1)/2;
        dv = (v2-v1)/2;

        for i = 1:N
            for j = 1:N

                dNu = dersbasisfuns( ...
                    iu,u_sample(i),Nurbs2D_stiff.order{k}(1)-1,2, ...
                    Nurbs2D_stiff.knots.U{k});

                dNv = dersbasisfuns( ...
                    iv,v_sample(j),Nurbs2D_stiff.order{k}(2)-1,2, ...
                    Nurbs2D_stiff.knots.V{k});

                [~,dS] = derRat2DBasisFuns( ...
                    dNu,dNv, ...
                    Nurbs2D_stiff.order{k}(1), ...
                    Nurbs2D_stiff.order{k}(2), ...
                    CP,2,2);

                nodeData(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);

                A1 = dS(:,2,1);
                A2 = dS(:,1,2);

                t1 = A1./norm(A1);

                A3 = cross(A1,A2);
                A3 = A3./norm(A3);

                t2 = cross(A3,t1);
                t2 = t2./norm(t2);

                % First fundamental form
                A = [A1,A2];
                g = A.'*A;

                if rcond(g) < 1E-12
                    error('semOpt2Dmesh_v1:DegenerateStiffenerMetric', ...
                        ['Degenerate stiffener surface metric at global element ' ...
                         '%d, sampling point %d.'],count_el,count);
                end

                % Contravariant basis
                Ac = A/g;

                % Signed second fundamental form in parametric coordinates
                A11 = dS(:,3,1);
                A12 = dS(:,2,2);
                A22 = dS(:,1,3);

                bParam = [dot(A11,A3), dot(A12,A3); ...
                          dot(A12,A3), dot(A22,A3)];

                % Curvature tensor in the local orthonormal basis {t1,t2}
                T = [t1,t2];
                C = g\(A.'*T);

                bLocal = C.'*bParam*C;
                bLocal = 0.5.*(bLocal+bLocal.');

                % Signed principal curvatures
                kappa = eig(bLocal);
                kappa = sort(real(kappa));

                JacMatData(:,:,count,count_el) = A;

                JacobianData(1,count,count_el) = ...
                    norm(cross(A1,A2))*du*dv;

                InvJacMatData(:,:,count,count_el) = ...
                    [dot(t1,Ac(:,1))/du, dot(t1,Ac(:,2))/dv; ...
                     dot(t2,Ac(:,1))/du, dot(t2,Ac(:,2))/dv];

                curvData(:,count,count_el) = kappa;
                curvTensorData(:,:,count,count_el) = bLocal;

                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';

                RData(:,:,count,count_el) = [t1,t2,A3];

                count = count + 1;
                count_node = count_node + 1;
            end
        end

        count_el = count_el + 1;
    end
end

%% Merge coincident plate-plate and plate-stiffener sampling points
TOL = 1E-5;

[nodes_sem,IA,IC] = uniquetol(nodeData,TOL,'ByRows',true);

elemNode = reshape(IC,N*N,nel).';

conn_sem = zeros(nel,shell_dof*N*N);

for d = 1:shell_dof
    conn_sem(:,d:shell_dof:end) = ...
        shell_dof*elemNode - (shell_dof-d);
end

%% Store SEM mesh and geometric quantities
semOpt2D.nodes = nodes_sem;
semOpt2D.conn = conn_sem;

semOpt2D.Jmat = JacMatData;
semOpt2D.J = JacobianData;
semOpt2D.InvJmat = InvJacMatData;

semOpt2D.Kappa = curvData;
semOpt2D.CurvTensor = curvTensorData;

semOpt2D.t1 = T1Data(IA,:);
semOpt2D.t2 = T2Data(IA,:);
semOpt2D.n  = NData(IA,:);

semOpt2D.R = RData;

%% Chebyshev operators: xi direction
space.a = -1;
space.b = 1;
space.N = N;

[FT_xi,BT_xi] = cheb(space);
D_xi = derivative(space);
V_xi = InnerProduct(space);

Q1_xi = BT_xi*D_xi*FT_xi;
Q2_xi = BT_xi*D_xi^2*FT_xi;

semOpt2D.FT_xi = FT_xi;
semOpt2D.BT_xi = BT_xi;
semOpt2D.D_xi = D_xi;
semOpt2D.V_xi = V_xi;
semOpt2D.Q1_xi = Q1_xi;
semOpt2D.Q2_xi = Q2_xi;

%% Chebyshev operators: eta direction
space.a = -1;
space.b = 1;
space.N = N;

[FT_eta,BT_eta] = cheb(space);
D_eta = derivative(space);
V_eta = InnerProduct(space);

Q1_eta = BT_eta*D_eta*FT_eta;
Q2_eta = BT_eta*D_eta^2*FT_eta;

semOpt2D.FT_eta = FT_eta;
semOpt2D.BT_eta = BT_eta;
semOpt2D.D_eta = D_eta;
semOpt2D.V_eta = V_eta;
semOpt2D.Q1_eta = Q1_eta;
semOpt2D.Q2_eta = Q2_eta;

%% Two-dimensional operators
semOpt2D.VD = kron(V_xi,V_eta);

semOpt2D.Q1xi = kron(Q1_xi,eye(N));
semOpt2D.Q1eta = kron(eye(N),Q1_eta);

semOpt2D.Q2xi = kron(Q2_xi,eye(N));
semOpt2D.Q2eta = kron(eye(N),Q2_eta);

semOpt2D.Qxieta = semOpt2D.Q2xi*semOpt2D.Q2eta;

end
