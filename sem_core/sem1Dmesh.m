function sem1D = sem1Dmesh(Nurbs1D,N,local_dof)
nel = 0; %Total number of elements
for k = 1:Nurbs1D.numpatch
    nel = nel + Nurbs1D.nel{k};
end
sem1D.nel = nel;
sem1D.N = N;
sem1D.local_dof = local_dof;
%
elData = zeros(N,3,nel);
nodeData = zeros(3*nel,3);
count_el = 1;
count_node = 1;
%
xi = lobat(N);
%
for k = 1:Nurbs1D.numpatch
    for el = 1:Nurbs1D.nel{k}
        iu = Nurbs1D.INC{k}(Nurbs1D.IEN{k}(1,el),1);
        u1 = Nurbs1D.knots.U{k}(iu);
        u2 = Nurbs1D.knots.U{k}(iu+1);
        u_sample = (0.5.*(1-xi).*u1+0.5.*(1+xi).*u2);
        count = 1;
        for i = 1:length(u_sample)
            dNu = dersbasisfuns(iu,u_sample(i),Nurbs1D.order{k}(1)-1,1,Nurbs1D.knots.U{k});
            CP = Nurbs1D.cPoints{k}(:,iu-Nurbs1D.order{k}(1)+1:iu);
            [~,dC] = derRat1DBasisFuns(dNu,Nurbs1D.order{k}(1),CP,1);
            elData(count,:,count_el) = dC(:,1)';
            nodeData(count_node,:) = dC(:,1)';
            count = count+1;
            count_node = count_node+1;
        end
        count_el = count_el+1;
    end
end
TOL = 0.00001; %---> Check!
nodes_sem = uniquetol(nodeData,TOL,'ByRows',true);
conn_sem = zeros(nel,local_dof*N);
for i = 1:nel
    for j = 1:N
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
sem1D.nodes = nodes_sem;
sem1D.conn = conn_sem;
sem1D.FT = FT_xi;
sem1D.BT = BT_xi;
sem1D.D = D_xi;
sem1D.V = V_xi;
sem1D.Q1 = Q1_xi;