function m_loc = localMass2D(sem2D,el)
%
ldof = sem2D.local_dof;
nconn = sem2D.conn(el,ldof:ldof:end)./ldof;
%
m_loc = zeros(ldof*sem2D.N*sem2D.N);
%
if sem2D.ET == 1
    %
    VD = sem2D.VD*diag(sem2D.J(nconn));
    %
    m_loc(1:ldof:end,1:ldof:end) = m_loc(1:ldof:end,1:ldof:end) + (sem2D.rho*sem2D.t).*VD;
    m_loc(2:ldof:end,2:ldof:end) = m_loc(2:ldof:end,2:ldof:end) + (sem2D.rho*sem2D.t^3/12).*VD;
    m_loc(3:ldof:end,3:ldof:end) = m_loc(3:ldof:end,3:ldof:end) + (sem2D.rho*sem2D.t^3/12).*VD;
end