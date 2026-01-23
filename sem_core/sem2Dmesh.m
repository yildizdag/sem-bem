function sem2D = sem2Dmesh(Nurbs2D,N,local_dof)
nel = 0;
for k = 1:Nurbs2D.numpatch
    nel = nel + Nurbs2D.nel{k};
end
sem2D.nel = nel;
sem2D.N = N;
sem2D.local_dof = local_dof;
%
nodeData = zeros(N*N*nel,3);
JacMatData = zeros(2,2,N*N*nel);
InvJacMatData = zeros(2,2,N*N*nel);
JacobianData = zeros(N*N*nel,1);
count_el = 1;
count_node = 1;
%
xi = lobat(N);
eta = lobat(N);
%
epsilon = 1E-6;
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
                dNu = dersbasisfuns(iu,u_sample(i),Nurbs2D.order{k}(1)-1,1,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(j),Nurbs2D.order{k}(2)-1,1,Nurbs2D.knots.V{k});
                % CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                [~,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{k}(1),Nurbs2D.order{k}(2),CP,1,1);
                nodeData(count_node,:) = epsilon.*(dS(:,1,1)'./epsilon);
% %                 elData(count,:,count_el) = epsilon.*(dS(:,1,1)'./epsilon);
                % % J2 = diag([Nurbs2D.knots.U{k}(iu+1)-Nurbs2D.knots.U{k}(iu),Nurbs2D.knots.V{k}(iv+1) - Nurbs2D.knots.V{k}(iv)])/2;
                % % JacMatData(:,:,count_node) = J2*[dS(1,2,1), dS(2,2,1); dS(1,1,2), dS(2,1,2)];
                % % InvJacMatData(:,:,count_node) = inv(J2*[dS(1,2,1), dS(2,2,1); dS(1,1,2), dS(2,1,2)]);
                % % JacobianData(count_node) = det(J2*[dS(1,2,1), dS(2,2,1); dS(1,1,2), dS(2,1,2)]);
% %                 du = (Nurbs2D.knots.U{k}(iu+1)-Nurbs2D.knots.U{k}(iu))/2;
% %                 dv = (Nurbs2D.knots.V{k}(iv+1)-Nurbs2D.knots.V{k}(iv))/2;
                %
                a = du * dS(1,2,1);
                b = du * dS(2,2,1);
                c = dv * dS(1,1,2);
                d = dv * dS(2,1,2);
                %
                detJ = a*d - b*c;
                %
                JacMatData(:,:,count_node) = [a b; c d];
                JacobianData(count_node)   = detJ;
                InvJacMatData(:,:,count_node) = (1/detJ).*[d -b; -c a];
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
% % TOL = 0.0001; %---> Check!
% % [~,IA] = uniquetol(nodeData,TOL,'ByRows',true);
% % IA = sort(IA); 
% % nodes_sem = nodeData(IA,:);
% % Jmat = JacMatData(:,:,IA);
% % InvJmat = InvJacMatData(:,:,IA);
% % J = JacobianData(IA);
% % conn_sem = zeros(nel,local_dof*N*N);
% % for i = 1:nel
% %     for j = 1:N*N
% %         node_id = find(ismembertol(nodes_sem, elData(j,:,i),TOL,'ByRows',true));
% %         for k = 1:local_dof
% %             conn_sem(i,local_dof*j-(local_dof-k)) = local_dof*node_id-(local_dof-k);
% %         end
% %     end
% % end
TOL = 1e-4; % keep your value for now
[nodes_sem, IA, IC] = uniquetol(nodeData, TOL, 'ByRows', true);
Jmat = JacMatData(:,:,IA);
InvJmat = InvJacMatData(:,:,IA);
J = JacobianData(IA);
elemNode = reshape(IC, N*N, nel).';
conn_sem = zeros(nel, local_dof*N*N);
for d = 1:local_dof
    conn_sem(:, d:local_dof:end) = local_dof*elemNode - (local_dof - d);
end
sem2D.nodes = nodes_sem;
sem2D.conn = conn_sem;
sem2D.Jmat = Jmat;
sem2D.J = J;
sem2D.InvJmat = InvJmat;
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