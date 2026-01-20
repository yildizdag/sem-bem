function k_loc = localStiffness1D(sem1D,el)
%
ldof = sem1D.local_dof;
nconn = sem1D.conn(el,ldof:ldof:end)./ldof;
if sem1D.ET == 1
    J = (sem1D.Q1)*sem1D.nodes(nconn,1);
    k_loc = (sem1D.EA).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
elseif sem1D.ET == 2
    k_loc = zeros(sem1D.local_dof*sem1D.N);
    J = (sem1D.Q1)*sem1D.nodes(nconn,1);
    k_loc(2:2:end,2:2:end) = k_loc(2:2:end,2:2:end) + (sem1D.EI).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(1:2:end,1:2:end) = k_loc(1:2:end,1:2:end) + (sem1D.kGA).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(1:2:end,2:2:end) = k_loc(1:2:end,2:2:end) - (sem1D.kGA).*((sem1D.Q1./J)'*sem1D.V).*J;
    k_loc(2:2:end,1:2:end) = k_loc(2:2:end,1:2:end) - (sem1D.kGA).*(sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(2:2:end,2:2:end) = k_loc(2:2:end,2:2:end) + (sem1D.kGA).*sem1D.V.*J;
end