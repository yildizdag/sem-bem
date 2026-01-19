function K = globalStiffness1D(sem1D)
storeSparse = (sem1D.N*sem1D.local_dof)*sem1D.nel;
%
I = zeros(storeSparse,1);
J = zeros(storeSparse,1);
K = zeros(storeSparse,1);
ntriplets = 0;
%
for el = 1:sem1D.nel
    el_conn = sem1D.conn(el,:);
    [~,i,s] = find(el_conn);
    k_loc = localStiffness1D(sem1D,el);
    for k = 1:numel(i)
        for l = 1:numel(i)
            ntriplets = ntriplets+1;
            I(ntriplets) = s(k);
            J(ntriplets) = s(l);
            K(ntriplets) = k_loc(i(k), i(l));
        end
    end
end
%
K = sparse(I(1:ntriplets),J(1:ntriplets),K(1:ntriplets),sem1D.dof,sem1D.dof);
K = (K+K')./2;