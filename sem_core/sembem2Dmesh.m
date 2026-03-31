function sembem2D = sembem2Dmesh(Nurbs2D,N,semPatch,bemPatch,shell_dof,fluid_dof)
%
numSEMpatch = size(semPatch,2);
numBEMpatch = size(bemPatch,2);
%
nel = 0;
for k = 1:numSEMpatch
    nel = nel + Nurbs2D.nel{k};
end
sembem2D.nel = nel;
sembem2D.N = N;
sembem2D.shell_dof = shell_dof;
sembem2D.fluid_dof = fluid_dof;
%
ntot = N*N*nel;
nodeData = zeros(ntot,3);
JacMatData = zeros(3,2,N*N,nel);
InvJacMatData = zeros(2,2,N*N,nel);
JacobianData = zeros(1,N*N,nel);
curvData = zeros(ntot,2);
%
T1Data = zeros(ntot,3);   % local tangent-1 (global components)
T2Data = zeros(ntot,3);   % local tangent-2 (global components)
NData  = zeros(ntot,3);   % local normal    (global components)
RData  = zeros(3,3,ntot); % optional: rotation matrix per sampling point
%
count_el = 1;
count_node = 1;
%
xi = lobat(N);
eta = lobat(N);
%
epsilon = 1E-6;
%
for k = 1:numSEMpatch
    for el = 1:Nurbs2D.nel{k}
        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);
        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);
        u_sample = (0.5.*(1-xi).*u1+0.5.*(1+xi).*u2);
        v_sample = (0.5.*(1-eta).*v1+0.5.*(1+eta).*v2);
        count = 1;
        CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
        du = (Nurbs2D.knots.U{k}(iu+1)-Nurbs2D.knots.U{k}(iu))/2;
        dv = (Nurbs2D.knots.V{k}(iv+1)-Nurbs2D.knots.V{k}(iv))/2;
        for i = 1:N
            for j = 1:N
                dNu = dersbasisfuns(iu,u_sample(i),Nurbs2D.order{k}(1)-1,2,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,2,Nurbs2D.knots.V{k});
                [~,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{k}(1),Nurbs2D.order{k}(2),CP,2,2);
                nodeData(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);
                %
                A1 = dS(:,2,1); A2 = dS(:,1,2);
                t1 = A1 ./ norm(A1);
                %
                A3 = cross(A1,A2)/norm(cross(A1,A2));
                %
                t2 = cross(A3,t1);
                t2 = t2 ./ norm(t2);
                %
                F1 = [A1 A2]'*[A1 A2];
                Ac = [A1, A2]/F1;
                F2 = [dot(dS(:,3,1),A3), dot(dS(:,2,2),A3)
                      dot(dS(:,2,2),A3), dot(dS(:,1,3),A3)];
                F = F1\F2;
                kappa = eig(F);
                %
                JacMatData(:,:,count,count_el) = [A1, A2];
                JacobianData(1,count,count_el) = norm(cross(A1,A2))*du*dv;
                InvJacMatData(:,:,count,count_el) = [dot(t1,Ac(:,1))/du dot(t1,Ac(:,2))/dv; dot(t2,Ac(:,1))/du dot(t2,Ac(:,2))/dv];
                curvData(count_node,:) = [abs(kappa(1)), abs(kappa(2))];
                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';
                RData(:,:,count_node) = [t1, t2, A3];
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
TOL = 1e-5;
[nodes_sem, IA, IC] = uniquetol(nodeData, TOL, 'ByRows', true);
Kappa = curvData(IA,:);
elemNode = reshape(IC, N*N, nel).';
conn_sem = zeros(nel, shell_dof*N*N);
for d = 1:shell_dof
    conn_sem(:, d:shell_dof:end) = shell_dof*elemNode - (shell_dof - d);
end
sembem2D.nodes = nodes_sem;
sembem2D.conn = conn_sem;
sembem2D.Jmat = JacMatData;
sembem2D.J = JacobianData;
sembem2D.InvJmat = InvJacMatData;
sembem2D.Kappa = Kappa;
sembem2D.t1 = T1Data(IA,:);
sembem2D.t2 = T2Data(IA,:);
sembem2D.n  = NData(IA,:);
sembem2D.R  = RData(:,:,IA); 
% xi-direction:
space.a=-1; space.b=1; space.N=N;
[FT_xi,BT_xi] = cheb(space);
D_xi = derivative(space);
V_xi = InnerProduct(space);
Q1_xi = BT_xi*D_xi*FT_xi;
Q2_xi = BT_xi*D_xi^2*FT_xi;
% Store:
sembem2D.FT_xi = FT_xi;
sembem2D.BT_xi = BT_xi;
sembem2D.D_xi = D_xi;
sembem2D.V_xi = V_xi;
sembem2D.Q1_xi = Q1_xi;
sembem2D.Q2_xi = Q2_xi;
% eta-direction:
space.a=-1; space.b=1; space.N=N;
[FT_eta,BT_eta] = cheb(space);
D_eta = derivative(space);
V_eta = InnerProduct(space);
Q1_eta = BT_eta*D_eta*FT_eta;
Q2_eta = BT_eta*D_eta^2*FT_eta;
% Store:
sembem2D.FT_eta = FT_eta;
sembem2D.BT_eta = BT_eta;
sembem2D.D_eta = D_eta;
sembem2D.V_eta = V_eta;
sembem2D.Q1_eta = Q1_eta;
sembem2D.Q2_eta = Q2_eta;
%
sembem2D.VD = kron(V_xi,V_eta);
sembem2D.Q1xi = kron(Q1_xi,eye(N));
sembem2D.Q1eta = kron(eye(N),Q1_eta);
%
sembem2D.Q2xi = kron(Q2_xi,eye(N));
sembem2D.Q2eta = kron(eye(N),Q2_eta);
sembem2D.Qxieta = sembem2D.Q2xi*sembem2D.Q2eta;
%
sembem2D.FT = Fxy_mapping(N,N,FT_xi,FT_eta);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
nelBEM = 0;
for k = numSEMpatch+1:numSEMpatch+numBEMpatch
    nelBEM = nelBEM + Nurbs2D.nel{k};
end
sembem2D.nelBEM = nelBEM;
%
ntotBEM = N*N*nelBEM;
nodeDataBEM = zeros(ntotBEM,3);
%
NDataBEM  = zeros(ntotBEM,3);
%
count_el = 1;
count_node = 1;
%
xi = linspace(-1,1,N);
eta = linspace(-1,1,N);
%
epsilon = 1E-6;
%
for k = numSEMpatch+1:numSEMpatch+numBEMpatch
    for el = 1:Nurbs2D.nel{k}
        iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
        iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
        u1 = Nurbs2D.knots.U{k}(iu);
        u2 = Nurbs2D.knots.U{k}(iu+1);
        v1 = Nurbs2D.knots.V{k}(iv);
        v2 = Nurbs2D.knots.V{k}(iv+1);
        u_sample = (0.5.*(1-xi).*u1+0.5.*(1+xi).*u2);
        v_sample = (0.5.*(1-eta).*v1+0.5.*(1+eta).*v2);
        count = 1;
        CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
        for i = 1:N
            for j = 1:N
                dNu = dersbasisfuns(iu,u_sample(i),Nurbs2D.order{k}(1)-1,2,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,2,Nurbs2D.knots.V{k});
                [~,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{k}(1),Nurbs2D.order{k}(2),CP,2,2);
                nodeDataBEM(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);
                %
                A1 = dS(:,2,1); A2 = dS(:,1,2);
                A3 = cross(A1,A2)/norm(cross(A1,A2));
                %
                NDataBEM(count_node,:)  = A3.';
                %
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
TOL = 1e-5;
[nodes_bem, IA, IC] = uniquetol(nodeDataBEM, TOL, 'ByRows', true);
elemNodeBEM = reshape(IC, N*N, nelBEM).';
conn_bem = zeros(nelBEM, fluid_dof*N*N);
for d = 1:fluid_dof
    conn_bem(:, d:fluid_dof:end) = fluid_dof*elemNodeBEM - (fluid_dof - d);
end
sembem2D.nodesBEM = nodes_bem;
sembem2D.connBEM = conn_bem;
sembem2D.nBEM = NDataBEM(IA,:);
