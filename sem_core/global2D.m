function [K,M] = global2D(sem2D)
%
el_dof = (sem2D.N*sem2D.N*sem2D.local_dof)^2;
storeSparse = el_dof*sem2D.nel;
%
SP = zeros(storeSparse,4);
%
ntriplets = 0;
%
for el = 1:sem2D.nel
    el_conn = sem2D.conn(el,:);
    [k_loc,m_loc] = local2D(sem2D,el);
    [I,J] = ndgrid(el_conn,el_conn);
    %
    ind = ntriplets + (1:el_dof);
    %
    SP(ind,:) = [I(:),J(:),k_loc(:),m_loc(:)];
    %
    ntriplets = ntriplets + (sem2D.N*sem2D.N*sem2D.local_dof)^2; 
end
%
K = sparse(SP(1:ntriplets,1),SP(1:ntriplets,2),SP(1:ntriplets,3),sem2D.dof,sem2D.dof);
M = sparse(SP(1:ntriplets,1),SP(1:ntriplets,2),SP(1:ntriplets,4),sem2D.dof,sem2D.dof);
%