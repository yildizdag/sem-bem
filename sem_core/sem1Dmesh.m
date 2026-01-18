function sem1D = sem1Dmesh(Nurbs1D,N)
nel = 0; %Total num of elements
for k = 1:Nurbs1D.numpatch
    nel = nel + Nurbs1D.nel{k};
end
sem1D.nel = nel;
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
        for i = 1:length(u_sample)
            dNu = dersbasisfuns(iu,u_sample(i),Nurbs1D.order{k}(1)-1,1,Nurbs1D.knots.U{k});
            CP = Nurbs1D.cPoints{k}(:,iu-Nurbs1D.order{k}(1)+1:iu);
            [dR,dC] = derRat1DBasisFuns(dNu,Nurbs1D.order{k}(1),CP,1);
        end
    end
end