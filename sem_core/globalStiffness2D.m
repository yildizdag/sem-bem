function K = globalStiffness2D(sem2D)
storeSparse = (sem2D.N*sem2D.N*sem2D.local_dof)^2*sem2D.nel;
%
I = zeros(storeSparse,1);
J = zeros(storeSparse,1);
K = zeros(storeSparse,1);
ntriplets = 0;
%
for el = 1:sem2D.nel
    el_conn = sem2D.conn(el,:);
    [~,i,s] = find(el_conn);
    k_loc = localStiffness2D(sem2D,el);
    for k = 1:numel(i)
        for l = 1:numel(i)
            ntriplets = ntriplets+1;
            I(ntriplets) = s(k);
            J(ntriplets) = s(l);
            K(ntriplets) = k_loc(i(k),i(l));
        end
    end
end
%
K = sparse(I(1:ntriplets),J(1:ntriplets),K(1:ntriplets),sem2D.dof,sem2D.dof);
K = (K+K')./2;