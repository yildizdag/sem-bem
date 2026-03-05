function [k_loc,m_loc] = local2D(sem2D,el)
%
n_dof = sem2D.shell_dof;
nconn = sem2D.conn(el,n_dof:n_dof:end)./n_dof;
n_el = sem2D.N*sem2D.N;
%
if sem2D.ET == 1 %-Flat Plate on x-y plane
    %
    k_loc = zeros(n_dof*n_el);
    m_loc = zeros(n_dof*n_el);
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
    %
elseif sem2D.ET == 2 %-Curved Shell (FSDT Shallow)
    %
    k_loc1 = zeros(5*n_el);
    m_loc1 = zeros(5*n_el);
    %
    VD = sem2D.VD * diag(sem2D.J(nconn));
    VD_beta1 = VD * diag(sem2D.Kappa(nconn,1));
    VD_beta2 = VD * diag(sem2D.Kappa(nconn,2));
    VD_beta1sq = VD * diag(sem2D.Kappa(nconn,1).^2);
    VD_beta2sq = VD * diag(sem2D.Kappa(nconn,2).^2);
    VD_beta1beta2 = VD * diag(sem2D.Kappa(nconn,1).*sem2D.Kappa(nconn,2));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,nconn),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(2,1,nconn),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(1,2,nconn),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,nconn),n_el,1).*sem2D.Q1eta;
    %
    QDx = QDxi_dxidx + QDeta_detadx;
    QDy = QDxi_dxidy + QDeta_detady;
    %
    C1 = (sem2D.t^3/12);
    C2 = sem2D.t;
    k_sc = 5/6;
    %
    k_loc1(1:5:end,1:5:end) = k_loc1(1:5:end,1:5:end) + sem2D.lame*C2*QDx'*VD*QDx;
    k_loc1(1:5:end,3:5:end) = k_loc1(1:5:end,3:5:end) + sem2D.lame*C2*QDx'*VD_beta1;
    k_loc1(1:5:end,2:5:end) = k_loc1(1:5:end,2:5:end) + sem2D.nu*sem2D.lame*C2*QDx'*VD*QDy;
    k_loc1(1:5:end,3:5:end) = k_loc1(1:5:end,3:5:end) + sem2D.nu*sem2D.lame*C2*QDx'*VD_beta2;
    k_loc1(1:5:end,1:5:end) = k_loc1(1:5:end,1:5:end) + sem2D.G*C2*QDy'*VD*QDy;
    k_loc1(1:5:end,2:5:end) = k_loc1(1:5:end,2:5:end) + sem2D.G*C2*QDy'*VD*QDx;
    k_loc1(1:5:end,1:5:end) = k_loc1(1:5:end,1:5:end) + k_sc*sem2D.G*C2*VD_beta1sq;
    k_loc1(1:5:end,3:5:end) = k_loc1(1:5:end,3:5:end) - k_sc*sem2D.G*C2*VD_beta1*QDx;
    k_loc1(1:5:end,4:5:end) = k_loc1(1:5:end,4:5:end) - k_sc*sem2D.G*C2*VD_beta1;
    %
    k_loc1(2:5:end,1:5:end) = k_loc1(2:5:end,1:5:end) + sem2D.nu*sem2D.lame*C2*QDy'*VD*QDx;
    k_loc1(2:5:end,3:5:end) = k_loc1(2:5:end,3:5:end) + sem2D.nu*sem2D.lame*C2*QDy'*VD_beta1;
    k_loc1(2:5:end,2:5:end) = k_loc1(2:5:end,2:5:end) + sem2D.lame*C2*QDy'*VD*QDy;
    k_loc1(2:5:end,3:5:end) = k_loc1(2:5:end,3:5:end) + sem2D.lame*C2*QDy'*VD_beta2;
    k_loc1(2:5:end,1:5:end) = k_loc1(2:5:end,1:5:end) + sem2D.G*C2*QDx'*VD*QDy;
    k_loc1(2:5:end,2:5:end) = k_loc1(2:5:end,2:5:end) + sem2D.G*C2*QDx'*VD*QDx;
    k_loc1(2:5:end,2:5:end) = k_loc1(2:5:end,2:5:end) + k_sc*sem2D.G*C2*VD_beta2sq;
    k_loc1(2:5:end,3:5:end) = k_loc1(2:5:end,3:5:end) - k_sc*sem2D.G*C2*VD_beta2*QDy;
    k_loc1(2:5:end,5:5:end) = k_loc1(2:5:end,5:5:end) - k_sc*sem2D.G*C2*VD_beta2;
    %
    k_loc1(3:5:end,1:5:end) = k_loc1(3:5:end,1:5:end) + sem2D.lame*C2*VD_beta1*QDx;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + sem2D.lame*C2*VD_beta1sq;
    k_loc1(3:5:end,2:5:end) = k_loc1(3:5:end,2:5:end) + sem2D.nu*sem2D.lame*C2*VD_beta1*QDy;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + sem2D.nu*sem2D.lame*C2*VD_beta1beta2;
    k_loc1(3:5:end,1:5:end) = k_loc1(3:5:end,1:5:end) + sem2D.nu*sem2D.lame*C2*VD_beta2*QDx;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + sem2D.nu*sem2D.lame*C2*VD_beta1beta2;
    k_loc1(3:5:end,2:5:end) = k_loc1(3:5:end,2:5:end) + sem2D.lame*C2*VD_beta2*QDy;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + sem2D.lame*C2*VD_beta2sq;
    k_loc1(3:5:end,1:5:end) = k_loc1(3:5:end,1:5:end) - k_sc*sem2D.G*C2*QDx'*VD_beta1;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + k_sc*sem2D.G*C2*QDx'*VD*QDx;
    k_loc1(3:5:end,4:5:end) = k_loc1(3:5:end,4:5:end) + k_sc*sem2D.G*C2*QDx'*VD;
    k_loc1(3:5:end,2:5:end) = k_loc1(3:5:end,2:5:end) - k_sc*sem2D.G*C2*QDy'*VD_beta2;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + k_sc*sem2D.G*C2*QDy'*VD*QDy;
    k_loc1(3:5:end,5:5:end) = k_loc1(3:5:end,5:5:end) + k_sc*sem2D.G*C2*QDy'*VD;
    %
    k_loc1(4:5:end,4:5:end) = k_loc1(4:5:end,4:5:end) + sem2D.lame*C1*QDx'*VD*QDx;
    k_loc1(4:5:end,5:5:end) = k_loc1(4:5:end,5:5:end) + sem2D.nu*sem2D.lame*C1*QDx'*VD*QDy;
    k_loc1(4:5:end,4:5:end) = k_loc1(4:5:end,4:5:end) + sem2D.G*C1*QDy'*VD*QDy;
    k_loc1(4:5:end,5:5:end) = k_loc1(4:5:end,5:5:end) + sem2D.G*C1*QDy'*VD*QDx;
    k_loc1(4:5:end,1:5:end) = k_loc1(4:5:end,1:5:end) - k_sc*sem2D.G*C2*VD_beta1;
    k_loc1(4:5:end,3:5:end) = k_loc1(4:5:end,3:5:end) + k_sc*sem2D.G*C2*VD*QDx;
    k_loc1(4:5:end,4:5:end) = k_loc1(4:5:end,4:5:end) + k_sc*sem2D.G*C2*VD;
    %
    k_loc1(5:5:end,4:5:end) = k_loc1(5:5:end,4:5:end) + sem2D.nu*sem2D.lame*C1*QDy'*VD*QDx;
    k_loc1(5:5:end,5:5:end) = k_loc1(5:5:end,5:5:end) + sem2D.lame*C1*QDy'*VD*QDy;
    k_loc1(5:5:end,4:5:end) = k_loc1(5:5:end,4:5:end) + sem2D.G*C1*QDx'*VD*QDy;
    k_loc1(5:5:end,5:5:end) = k_loc1(5:5:end,5:5:end) + sem2D.G*C1*QDx'*VD*QDx;
    k_loc1(5:5:end,2:5:end) = k_loc1(5:5:end,2:5:end) - k_sc*sem2D.G*C2*VD_beta2;
    k_loc1(5:5:end,3:5:end) = k_loc1(5:5:end,3:5:end) + k_sc*sem2D.G*C2*VD*QDy;
    k_loc1(5:5:end,5:5:end) = k_loc1(5:5:end,5:5:end) + k_sc*sem2D.G*C2*VD;
    %
    m_loc1(1:5:end,1:5:end) = m_loc1(1:5:end,1:5:end) + sem2D.rho*C2*VD;
    m_loc1(2:5:end,2:5:end) = m_loc1(2:5:end,2:5:end) + sem2D.rho*C2*VD;
    m_loc1(3:5:end,3:5:end) = m_loc1(3:5:end,3:5:end) + sem2D.rho*C2*VD;
    m_loc1(4:5:end,4:5:end) = m_loc1(4:5:end,4:5:end) + sem2D.rho*C1*VD;
    m_loc1(5:5:end,5:5:end) = m_loc1(5:5:end,5:5:end) + sem2D.rho*C1*VD;
    %
    % indices of the 5 dofs inside the 6-dof layout, node-wise
    idx5 = zeros(5*n_el,1);
    for a = 1:n_el
        idx5( (a-1)*5 + (1:5) ) = (a-1)*6 + (1:5);   % place [u v w thx thy] into slots 1..5
    end
    idx6 = 6:6:(6*n_el);
    %
    minK = min(abs(diag(k_loc1)));
    minM = min(abs(diag(m_loc1)));
    %
    epsK = 1e-8 * minK;
    epsM = 1e-8 * minM;
    %
    k_loc = zeros(6*n_el);
    m_loc = zeros(6*n_el);
    k_loc(idx5,idx5) = k_loc1;
    m_loc(idx5,idx5) = m_loc1;
    k_loc(idx6,idx6) = epsK * eye(n_el);
    m_loc(idx6,idx6) = epsM * eye(n_el);
    %
    T = zeros(6*n_el);
    for a = 1:n_el
        R = sem2D.R(:,:, nconn(a));
        Tnode = blkdiag(R, R);
        ia = (a-1)*6 + (1:6);
        T(ia, ia) = Tnode;
    end
    k_loc = transpose(T) * k_loc * T;
    m_loc = transpose(T) * m_loc * T;
end