function sem2D = sem2Dmesh(Nurbs2D,N,local_dof)
nel = 0;
for k = 1:Nurbs2D.numpatch
    nel = nel + Nurbs2D.nel{k};
end
sem2D.nel = nel;
sem2D.N = N;
sem2D.local_dof = local_dof;
%
elData = zeros(N*N,3,nel);
nodeData = zeros(N*N*nel,3);
count_el = 1;
count_node = 1;
%
xi = lobat(N);
eta = lobat(N);
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
        for i = 1:length(v_sample)
            for j = 1:length(u_sample)
                dNu = dersbasisfuns(iu,u_sample(j),Nurbs2D.order{k}(1)-1,1,Nurbs2D.knots.U{k});
                dNv = dersbasisfuns(iv,v_sample(i),Nurbs2D.order{k}(2)-1,1,Nurbs2D.knots.V{k});
                CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                [~,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{k}(1),Nurbs2D.order{k}(2),CP,1,1);
                elData(count,:,count_el) = dS(:,1)';
                nodeData(count_node,:) = dS(:,1)';
                count = count+1;
                count_node = count_node+1;
            end
        end
        count_el = count_el+1;
    end
end
TOL = 0.00001; %---> Check!
[~,IA] = uniquetol(nodeData,TOL,'ByRows',true);
IA = sort(IA); 
nodes_sem = nodeData(IA,:);
conn_sem = zeros(nel,local_dof*N*N);
for i = 1:nel
    for j = 1:N*N
        node_id = find(ismembertol(nodes_sem, elData(j,:,i),TOL,'ByRows',true));
        for k = 1:local_dof
            conn_sem(i,local_dof*j-(local_dof-k)) = local_dof*node_id-(local_dof-k);
        end
    end
end
%
space.a=-1; space.b=1; space.N=N;
%
[FT_xi,BT_xi] = cheb(space);
D_xi = derivative(space);
V_xi = InnerProduct(space);
Q1_xi = BT_xi*D_xi*FT_xi;
%
sem2D.nodes = nodes_sem;
sem2D.conn = conn_sem;
sem2D.FT_xi = FT_xi;
sem2D.BT_xi = BT_xi;
sem2D.D_xi = D_xi;
sem2D.V_xi = V_xi;
sem2D.Q1_xi = Q1_xi;
%
space.a=-1; space.b=1; space.N=N;
%
[FT_eta,BT_eta] = cheb(space);
D_eta = derivative(space);
V_eta = InnerProduct(space);
Q1_eta = BT_eta*D_eta*FT_eta;
%
sem2D.nodes = nodes_sem;
sem2D.conn = conn_sem;
sem2D.FT_eta = FT_eta;
sem2D.BT_eta = BT_eta;
sem2D.D_eta = D_eta;
sem2D.V_eta = V_eta;
sem2D.Q1_eta = Q1_eta;