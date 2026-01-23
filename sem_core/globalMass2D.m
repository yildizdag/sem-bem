function M = globalMass2D(sem2D)
storeSparse = (sem2D.N*sem2D.N*sem2D.local_dof)^2*sem2D.nel;
%
P = zeros(storeSparse,1);
R = zeros(storeSparse,1);
S = zeros(storeSparse,1);
ntriplets = 0;
%
for el = 1:sem2D.nel
    el_conn = sem2D.conn(el,:);
    [~,i,s] = find(el_conn);
    m_loc = localMass2D(sem2D,el);
    for k = 1:numel(i)
        for l = 1:numel(i)
            ntriplets = ntriplets+1;
            P(ntriplets) = s(k);
            R(ntriplets) = s(l);
            S(ntriplets) = m_loc(i(k),i(l));
        end
    end
end
%
M = sparse(P(1:ntriplets),R(1:ntriplets),S(1:ntriplets),sem2D.dof,sem2D.dof);
%M = (M+M')./2;