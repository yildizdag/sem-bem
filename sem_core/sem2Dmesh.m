function sem2D = sem2Dmesh(Nurbs2D,N,shell_dof)
nel = 0;
for k = 1:Nurbs2D.numpatch
    nel = nel + Nurbs2D.nel{k};
end
sem2D.nel = nel;
sem2D.N = N;
sem2D.shell_dof = shell_dof;
%
ntot = N*N*nel;
nodeData = zeros(ntot,3);
JacMatData = zeros(3,2,ntot);
InvJacMatData = zeros(2,2,ntot);
JacobianData = zeros(ntot,1);
curvData = zeros(ntot,2);
%
T1Data = zeros(ntot,3);   % local tangent-1 (global components)
T2Data = zeros(ntot,3);   % local tangent-2 (global components)
NData  = zeros(ntot,3);   % local normal   (global components)
RData  = zeros(3,3,ntot); % optional: rotation matrix per sampling point
%
count_el = 1;
count_node = 1;
%
xi = lobat(N);
eta = lobat(N);
%
epsilon = 1E-4;
%
for k = 1:Nurbs2D.numpatch
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
                t2 = A2 ./ norm(A2);
                A3 = cross(A1,A2)/norm(cross(A1,A2));
                F1 = [A1 A2]'*[A1 A2];
                Ac = [A1, A2]/F1;
                F2 = [dot(dS(:,3,1),A3), dot(dS(:,2,2),A3)
                      dot(dS(:,2,2),A3), dot(dS(:,1,3),A3)];
                F = F1\F2;
                kappa = eig(F);
                %
                %a = du * dS(1,2,1);
                %b = du * dS(2,2,1);
                %c = dv * dS(1,1,2);
                %d = dv * dS(2,1,2);
                %
                %detJ = a*d - b*c;
                %
                JacMatData(:,:,count_node) = [A1, A2];
                JacobianData(count_node)   = norm(cross(A1,A2))*du*dv;
                InvJacMatData(:,:,count_node) = [dot(t1,Ac(:,1))/du dot(t1,Ac(:,2))/dv; dot(t2,Ac(:,1))/du dot(t2,Ac(:,2))/dv];
                curvData(count_node,:) = [abs(kappa(1)), abs(kappa(2))];
                T1Data(count_node,:) = t1.';
                T2Data(count_node,:) = t2.';
                NData(count_node,:)  = A3.';
                RData(:,:,count_node) = [t1 t2 A3];
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
TOL = 1e-4;
[nodes_sem, IA, IC] = uniquetol(nodeData, TOL, 'ByRows', true);
Jmat = JacMatData(:,:,IA);
InvJmat = InvJacMatData(:,:,IA);
J = JacobianData(IA);
Kappa = curvData(IA,:);
elemNode = reshape(IC, N*N, nel).';
conn_sem = zeros(nel, shell_dof*N*N);
for d = 1:shell_dof
    conn_sem(:, d:shell_dof:end) = shell_dof*elemNode - (shell_dof - d);
end
sem2D.nodes = nodes_sem;
sem2D.conn = conn_sem;
sem2D.Jmat = Jmat;
sem2D.J = J;
sem2D.InvJmat = InvJmat;
sem2D.Kappa = Kappa;
sem2D.t1 = T1Data(IA,:);
sem2D.t2 = T2Data(IA,:);
sem2D.n  = NData(IA,:);
sem2D.R  = RData(:,:,IA); 
% xi-direction:
space.a=-1; space.b=1; space.N=N;
[FT_xi,BT_xi] = cheb(space);
D_xi = derivative(space);
V_xi = InnerProduct(space);
Q1_xi = BT_xi*D_xi*FT_xi;
%
sem2D.FT_xi = FT_xi;
sem2D.BT_xi = BT_xi;
sem2D.D_xi = D_xi;
sem2D.V_xi = V_xi;
sem2D.Q1_xi = Q1_xi;
% eta-direction:
space.a=-1; space.b=1; space.N=N;
%
[FT_eta,BT_eta] = cheb(space);
D_eta = derivative(space);
V_eta = InnerProduct(space);
Q1_eta = BT_eta*D_eta*FT_eta;
%
sem2D.FT_eta = FT_eta;
sem2D.BT_eta = BT_eta;
sem2D.D_eta = D_eta;
sem2D.V_eta = V_eta;
sem2D.Q1_eta = Q1_eta;
%
sem2D.VD = kron(V_xi,V_eta);
sem2D.Q1xi = kron(Q1_xi,eye(N));
sem2D.Q1eta = kron(eye(N),Q1_eta);