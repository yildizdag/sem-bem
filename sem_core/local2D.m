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
    VD = sem2D.VD * diag(sem2D.J(1,nconn,el));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,nconn,el),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(2,1,nconn,el),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(1,2,nconn,el),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,nconn,el),n_el,1).*sem2D.Q1eta;
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
    epsK = 1e-4 * minK;
    epsM = 1e-4 * minM;
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
        Tnode = blkdiag(R.', R.');
        ia = (a-1)*6 + (1:6);
        T(ia, ia) = Tnode;
    end
    k_loc = transpose(T) * k_loc * T;
    m_loc = transpose(T) * m_loc * T;
    %
elseif sem2D.ET == 3
    %
    rho = sem2D.matProp(1);
    E11 = sem2D.matProp(2);
    E22 = sem2D.matProp(3);
    %
    G12 = sem2D.matProp(4);
    G31 = sem2D.matProp(5);
    G23 = sem2D.matProp(6);
    %
    nu12 = sem2D.matProp(7);
    nu21 = sem2D.matProp(8);
    %
    Q11 = E11/(1-nu12*nu21);
    Q22 = E22/(1-nu12*nu21);
    Q12 = nu12*E22/(1-nu12*nu21);
    Q66 = G12;
    Q44 = G23;
    Q55 = G31;
    %-Material Invariants
    matinv.U1 = 1/8*(3*Q11+3*Q22+2*Q12+4*Q66);
    matinv.U2 = 1/2*(Q11-Q22);
    matinv.U3 = 1/8*(Q11+Q22-2*Q12-4*Q66);
    matinv.U4 = 1/8*(Q11+Q22+6*Q12-4*Q66);
    matinv.U5 = 1/8*(Q11+Q22-2*Q12+4*Q66);

    %material invariant for shear partDS
    matinv.U11 = 1/2 * (Q44 + Q55);
    matinv.U22 = 1/2 * (Q44 - Q55);

    [A, B, D, G] = Stacking_LP_ABD(sem2D.stSeq,matinv,sem2D.t);
    %
    A11 = A(1,1); A22 = A(2,1); A12 = A(3,1);
    A66 = A(4,1); A16 = A(5,1); A26 = A(6,1);
    %
    A55 = G(1,1); A44 = G(2,1); A45 = G(3,1);
    %
    B11 = B(1,1); B22 = B(2,1); B12 = B(3,1);
    B66 = B(4,1); B16 = B(5,1); B26 = B(6,1);
    %
    D11 = D(1,1); D22 = D(2,1); D12 = D(3,1);
    D66 = D(4,1); D16 = D(5,1); D26 = D(6,1);
    %
    Kc = 5/6;
    %
    k_loc = zeros(5*n_el);
    m_loc = zeros(5*n_el);
    %
    VD = sem2D.VD * diag(sem2D.J(nconn));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,nconn),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(1,2,nconn),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(1,2,nconn),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,nconn),n_el,1).*sem2D.Q1eta;
    %
    QDx = QDxi_dxidx + QDeta_detadx;
    QDy = QDxi_dxidy + QDeta_detady;
    %
    k_loc(1:5:end,1:5:end) = k_loc(1:5:end,1:5:end) + A11*QDx'*VD*QDx + A16*QDy'*VD*QDx + A16*QDx'*VD*QDy + A66*QDy'*VD*QDy;
    k_loc(1:5:end,2:5:end) = k_loc(1:5:end,2:5:end) + A12*QDx'*VD*QDy + A26*QDy'*VD*QDy + A16*QDx'*VD*QDx + A66*QDy'*VD*QDx;
    %k_loc(1:5:end,3:5:end) = k_loc(1:5:end,3:5:end) + 0;
    k_loc(1:5:end,4:5:end) = k_loc(1:5:end,4:5:end) + B11*QDx'*VD*QDx + B16*QDy'*VD*QDx + B16*QDx'*VD*QDy + B66*QDy'*VD*QDy;
    k_loc(1:5:end,5:5:end) = k_loc(1:5:end,5:5:end) + B12*QDx'*VD*QDy + B26*QDy'*VD*QDy + B16*QDx'*VD*QDx + B66*QDy'*VD*QDx;
    %
    k_loc(2:5:end,1:5:end) = k_loc(2:5:end,1:5:end) + A12*QDy'*VD*QDx + A16*QDx'*VD*QDx + A26*QDy'*VD*QDy + A66*QDx'*VD*QDy;
    k_loc(2:5:end,2:5:end) = k_loc(2:5:end,2:5:end) + A22*QDy'*VD*QDy + A26*QDx'*VD*QDy + A26*QDy'*VD*QDx + A66*QDx'*VD*QDx;
    %k_loc(2:5:end,3:5:end) = 0;
    k_loc(2:5:end,4:5:end) = k_loc(2:5:end,4:5:end) + B12*QDy'*VD*QDx + B16*QDx'*VD*QDx + B26*QDy'*VD*QDy + B66*QDx'*VD*QDy;
    k_loc(2:5:end,5:5:end) = k_loc(2:5:end,5:5:end) + B22*QDy'*VD*QDy + B26*QDx'*VD*QDy + B26*QDy'*VD*QDx + B66*QDx'*VD*QDx;
    %
    %k_loc(3:5:end,1:5:end) = 0;
    %k_loc(3:5:end,2:5:end) = 0;
    k_loc(3:5:end,3:5:end) = k_loc(3:5:end,3:5:end) + Kc*A55*QDx'*VD*QDx + Kc*A45*QDy'*VD*QDx + Kc*A45*QDx'*VD*QDy + Kc*A44*QDy'*VD*QDy;
    k_loc(3:5:end,4:5:end) = k_loc(3:5:end,4:5:end) + Kc*A55*QDx'*VD + Kc*A45*QDy'*VD;
    k_loc(3:5:end,5:5:end) = k_loc(3:5:end,5:5:end) + Kc*A45*QDx'*VD + Kc*A44*QDy'*VD;
    %
    k_loc(4:5:end,1:5:end) = k_loc(4:5:end,1:5:end) + B11*QDx'*VD*QDx + B16*QDy'*VD*QDx + B16*QDx'*VD*QDy + B66*QDy'*VD*QDy ;
    k_loc(4:5:end,2:5:end) = k_loc(4:5:end,2:5:end) + B12*QDx'*VD*QDy + B26*QDy'*VD*QDy + B16*QDx'*VD*QDx + B66*QDy'*VD*QDx;
    k_loc(4:5:end,3:5:end) = k_loc(4:5:end,3:5:end) + Kc*A55*VD*QDx + Kc*A45*VD*QDy;
    k_loc(4:5:end,4:5:end) = k_loc(4:5:end,4:5:end) + D11*QDx'*VD*QDx + D16*QDy'*VD*QDx + D16*QDx'*VD*QDy + D66*QDy'*VD*QDy + Kc*A55*VD;
    k_loc(4:5:end,5:5:end) = k_loc(4:5:end,5:5:end) + D12*QDx'*VD*QDy + D26*QDy'*VD*QDy + D16*QDx'*VD*QDx + D66*QDy'*VD*QDx + Kc*A45*VD;
    %
    k_loc(5:5:end,1:5:end) = k_loc(5:5:end,1:5:end) + B12*QDy'*VD*QDx + B16*QDx'*VD*QDx + B26*QDy'*VD*QDy + B66*QDx'*VD*QDy;
    k_loc(5:5:end,2:5:end) = k_loc(5:5:end,2:5:end) + B22*QDy'*VD*QDy + B26*QDx'*VD*QDy + B26*QDy'*VD*QDx + B66*QDx'*VD*QDx;
    k_loc(5:5:end,3:5:end) = k_loc(5:5:end,3:5:end) + Kc*A45*VD*QDx + Kc*A44*VD*QDy;
    k_loc(5:5:end,4:5:end) = k_loc(5:5:end,4:5:end) + D12*QDy'*VD*QDx + D16*QDx'*VD*QDx + D26*QDy'*VD*QDy + D66*QDx'*VD*QDy + Kc*A45*VD;
    k_loc(5:5:end,5:5:end) = k_loc(5:5:end,5:5:end) + D22*QDy'*VD*QDy + D26*QDx'*VD*QDy + D26*QDy'*VD*QDx + D66*QDx'*VD*QDx + Kc*A44*VD;
    %
    m_loc(1:5:end,1:5:end) = m_loc(1:5:end,1:5:end) + rho*sem2D.t*VD;
    m_loc(2:5:end,2:5:end) = m_loc(2:5:end,2:5:end) + rho*sem2D.t*VD;
    m_loc(3:5:end,3:5:end) = m_loc(3:5:end,3:5:end) + rho*sem2D.t*VD;
    %
    m_loc(4:5:end,4:5:end) = m_loc(4:5:end,4:5:end) + (rho*sem2D.t^3/12)*VD;
    m_loc(5:5:end,5:5:end) = m_loc(5:5:end,5:5:end) + (rho*sem2D.t^3/12)*VD;

    k_loc = 0.5.*(k_loc + transpose(k_loc));
    m_loc = 0.5.*(m_loc + transpose(m_loc));
end