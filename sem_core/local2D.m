function [k_loc,m_loc] = local2D(sem2D,el)
%
n_dof = sem2D.local_dof;
nconn = sem2D.conn(el,n_dof:n_dof:end)./n_dof;
n_el = sem2D.N*sem2D.N;
%
k_loc = zeros(n_dof*n_el);
m_loc = zeros(n_dof*n_el);
%
if sem2D.ET == 1
    %
    VD = sem2D.VD * diag(sem2D.J(nconn));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,nconn),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(2,1,nconn),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(1,2,nconn),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,nconn),n_el,1).*sem2D.Q1eta;
    %
    QDx = QDxi_dxidx + QDeta_detadx;
    QDy = QDxi_dxidy + QDeta_detady;
    %
    k_loc(2:n_dof:end,2:n_dof:end) = k_loc(2:n_dof:end,2:n_dof:end) + sem2D.D.*(transpose(QDx)*VD*QDx) + sem2D.Ds.*VD + (sem2D.D*(1-sem2D.nu)/2).*(transpose(QDy)*VD*QDy);
    k_loc(2:n_dof:end,3:n_dof:end) = k_loc(2:n_dof:end,3:n_dof:end) + (sem2D.D*sem2D.nu).*(transpose(QDx)*VD*QDy) + (sem2D.D*(1-sem2D.nu)/2).*(transpose(QDy)*VD*QDx);
    k_loc(3:n_dof:end,2:n_dof:end) = k_loc(3:n_dof:end,2:n_dof:end) + (sem2D.D*sem2D.nu).*(transpose(QDy)*VD*QDx) + (sem2D.D*(1-sem2D.nu)/2).*(transpose(QDx)*VD*QDy);
    k_loc(3:n_dof:end,3:n_dof:end) = k_loc(3:n_dof:end,3:n_dof:end) + sem2D.D.*(transpose(QDy)*VD*QDy) + sem2D.Ds.*VD + (sem2D.D*(1-sem2D.nu)/2).*(transpose(QDx)*VD*QDx);
    k_loc(2:n_dof:end,1:n_dof:end) = k_loc(2:n_dof:end,1:n_dof:end) - sem2D.Ds.*(VD*QDx);
    k_loc(1:n_dof:end,2:n_dof:end) = k_loc(1:n_dof:end,2:n_dof:end) - sem2D.Ds.*(transpose(QDx)*VD);
    k_loc(1:n_dof:end,1:n_dof:end) = k_loc(1:n_dof:end,1:n_dof:end) + sem2D.Ds.*(transpose(QDx)*VD*QDx) + sem2D.Ds.*(transpose(QDy)*VD*QDy);
    k_loc(3:n_dof:end,1:n_dof:end) = k_loc(3:n_dof:end,1:n_dof:end) - sem2D.Ds.*(VD*QDy);
    k_loc(1:n_dof:end,3:n_dof:end) = k_loc(1:n_dof:end,3:n_dof:end) - sem2D.Ds.*(transpose(QDy)*VD);
    %
    m_loc(1:n_dof:end,1:n_dof:end) = m_loc(1:n_dof:end,1:n_dof:end) + (sem2D.rho*sem2D.t).*VD;
    m_loc(2:n_dof:end,2:n_dof:end) = m_loc(2:n_dof:end,2:n_dof:end) + (sem2D.rho*sem2D.t^3/12).*VD;
    m_loc(3:n_dof:end,3:n_dof:end) = m_loc(3:n_dof:end,3:n_dof:end) + (sem2D.rho*sem2D.t^3/12).*VD;
end