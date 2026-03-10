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
    A1  = sem2D.Q1xi * sem2D.nodes(nconn,:);
    A2  = sem2D.Q1eta * sem2D.nodes(nconn,:);
    %
    F1_11 = sum(A1.*A1,2);
    F1_12 = sum(A1.*A2,2);
    F1_22 = sum(A2.*A2,2);
    %
    A11   = sem2D.Q2xi * sem2D.nodes(nconn,:);
    A22   = sem2D.Q2eta * sem2D.nodes(nconn,:);
    A12   = sem2D.Qxieta * sem2D.nodes(nconn,:);
    %
    A3 = cross(A1,A2,2);
    A3 = A3 ./ vecnorm(A3,2,2);
    J = vecnorm(cross(A1,A2,2),2,2);
    %
    F2_11 = sum(A3.*A11,2);
    F2_12 = sum(A3.*A12,2);
    F2_22 = sum(A3.*A22,2);
    %
    Kappa  = zeros(n_el,2);
    t1     = zeros(n_el,3);
    t2     = zeros(n_el,3);
    %
    for a = 1:n_el
        %
        F1 = [F1_11(a) F1_12(a);
            F1_12(a) F1_22(a)];
        %
        F2 = [F2_11(a) F2_12(a);
            F2_12(a) F2_22(a)];
        %
        [V,D] = eig(F2,F1);
        %
        k = real(diag(D));
        [k,ord] = sort(k);
        %
        V = V(:,ord);
        % principal directions
        p1 = V(1,1)*A1(a,:) + V(2,1)*A2(a,:);
        p2 = V(1,2)*A1(a,:) + V(2,2)*A2(a,:);

        p1 = p1 / norm(p1);
        p2 = p2 / norm(p2);

        % enforce orthogonality
        n  = A3(a,:);
        p2 = cross(n,p1);
        p2 = p2 / norm(p2);

        Kappa(a,:) = k';

        t1(a,:) = p1;
        t2(a,:) = p2;
    end
    %
    VD = sem2D.VD * diag(J);
%     VD = sem2D.VD * diag(sem2D.J(nconn));
    VD_beta1 = VD * diag(Kappa(:,1));
    VD_beta2 = VD * diag(Kappa(:,2));
    VD_beta1sq = VD * diag(Kappa(:,1).^2);
    VD_beta2sq = VD * diag(Kappa(:,2).^2);
    VD_beta1beta2 = VD * diag(Kappa(:,1).*Kappa(:,2));
%     VD_beta1 = VD * diag(sem2D.Kappa(nconn,1));
%     VD_beta2 = VD * diag(sem2D.Kappa(nconn,2));
%     VD_beta1sq = VD * diag(sem2D.Kappa(nconn,1).^2);
%     VD_beta2sq = VD * diag(sem2D.Kappa(nconn,2).^2);
%     VD_beta1beta2 = VD * diag(sem2D.Kappa(nconn,1).*sem2D.Kappa(nconn,2));
    %
    %
    detF1 = F1_11.*F1_22 - F1_12.^2;
    %
    Ac1 = ( A1.*F1_22 - A2.*F1_12 ) ./ detF1;
    Ac2 = ( A2.*F1_11 - A1.*F1_12 ) ./ detF1;
    %
    InvJ11 = sum(t1.*Ac1,2);
    InvJ12 = sum(t1.*Ac2,2);
    InvJ21 = sum(t2.*Ac1,2);
    InvJ22 = sum(t2.*Ac2,2);
    %
%     QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,nconn),n_el,1).*sem2D.Q1xi;
%     QDxi_dxidy    = reshape(sem2D.InvJmat(2,1,nconn),n_el,1).*sem2D.Q1xi;
%     QDeta_detadx  = reshape(sem2D.InvJmat(1,2,nconn),n_el,1).*sem2D.Q1eta;
%     QDeta_detady  = reshape(sem2D.InvJmat(2,2,nconn),n_el,1).*sem2D.Q1eta;
    %
    R = zeros(3,3,n_el);
    %
    for a = 1:n_el
        R(:,:,a) = [t1(a,:); t2(a,:); A3(a,:)]';
    end
    %
    QDxi_dxidx    = reshape(InvJ11,n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(InvJ12,n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(InvJ21,n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(InvJ22,n_el,1).*sem2D.Q1eta;
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

end