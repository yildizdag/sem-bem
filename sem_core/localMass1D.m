function m_loc = localMass1D(sem1D,el)
%
ldof = sem1D.local_dof;
nconn = sem1D.conn(el,ldof:ldof:end)./ldof;
if sem1D.ET == 1
    J = (sem1D.Q1)*sem1D.nodes(nconn,1);
    m_loc = (sem1D.EA).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
elseif sem1D.ET == 2
    m_loc = zeros(sem1D.local_dof*sem1D.N);
    J = (sem1D.Q1)*sem1D.nodes(nconn,1);
    m_loc(2:2:end,2:2:end) = m_loc(2:2:end,2:2:end) + (sem1D.rho*sem1D.I).*(sem1D.V.*J);
    m_loc(1:2:end,1:2:end) = m_loc(1:2:end,1:2:end) + (sem1D.rho*sem1D.A).*(sem1D.V.*J);
end