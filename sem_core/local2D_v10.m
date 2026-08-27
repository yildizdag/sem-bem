function [k_loc,m_loc] = local2D_v10(sem2D,el)
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
    QDxi_dxidy    = reshape(sem2D.InvJmat(1,2,nconn,el),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(2,1,nconn,el),n_el,1).*sem2D.Q1eta;
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
elseif sem2D.ET == 2 %-Unified Reissner-Mindlin shell/stiffener formulation
    %
    % Global nodal DOFs at each SEM point:
    %   [ux uy uz rx ry rz]
    %
    % The formulation is orientation independent. Displacements and physical
    % rotation vectors are interpolated and differentiated in global Cartesian
    % coordinates. Their derivatives are then projected onto the pointwise
    % orthonormal surface frame [t1 t2 n]. Consequently, changes of the local
    % basis and the complete signed surface curvature (including b12) enter
    % through the geometry itself; no plate/stiffener switch is required.
    %
    if n_dof ~= 6
        error('local2D_v10:InvalidDOF', ...
            'The unified shell formulation requires six DOFs per SEM point.');
    end
    if ~isfield(sem2D,'CurvTensor')
        error('local2D_v10:MissingCurvatureTensor', ...
            ['sem2D.CurvTensor is missing. Generate the mesh with ', ...
             'sem2Dmesh_v1.m or a later compatible mesh routine.']);
    end
    %
    % Integration and physical derivative operators along the pointwise
    % orthonormal tangent directions t1 and t2.
    VD = sem2D.VD * diag(sem2D.J(1,:,el));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,:,el),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(1,2,:,el),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(2,1,:,el),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,:,el),n_el,1).*sem2D.Q1eta;
    %
    QDx = QDxi_dxidx + QDeta_detadx;
    QDy = QDxi_dxidy + QDeta_detady;
    %
    C1 = sem2D.t^3/12;
    C2 = sem2D.t;
    k_sc = 5/6;
    %
    % Extraction matrices for global translations and physical rotations.
    Eu = zeros(3*n_el,6*n_el);
    Er = zeros(3*n_el,6*n_el);
    Rblk = zeros(3*n_el);
    BetaOp = zeros(3*n_el);
    NormalOp = zeros(n_el,3*n_el);
    %
    for a = 1:n_el
        i3 = (a-1)*3 + (1:3);
        i6 = (a-1)*6 + (1:6);
        %
        Eu(i3,i6(1:3)) = eye(3);
        Er(i3,i6(4:6)) = eye(3);
        %
        Ra = sem2D.R(:,:,a,el);
        t1 = Ra(:,1);
        t2 = Ra(:,2);
        nv = Ra(:,3);
        % Re-orthonormalize defensively against small numerical drift.
        t1 = t1/norm(t1);
        nv = nv - t1*(t1.'*nv);
        nv = nv/norm(nv);
        t2 = cross(nv,t1);
        t2 = t2/norm(t2);
        Ra = [t1,t2,nv];
        Rblk(i3,i3) = Ra.';
        %
        % beta = theta x n = -skew(n)*theta.
        Sn = [      0, -nv(3),  nv(2); ...
               nv(3),       0, -nv(1); ...
              -nv(2),  nv(1),       0];
        BetaOp(i3,i3) = -Sn;
        NormalOp(a,i3) = nv.';
    end
    %
    Dx3 = kron(QDx,eye(3));
    Dy3 = kron(QDy,eye(3));
    %
    % Global displacement gradients, subsequently projected to [t1 t2 n].
    % Because differentiation is performed before projection, the connection
    % and signed-curvature terms associated with the varying frame are retained
    % implicitly and consistently, including the effect of b12.
    Ux = Rblk*Dx3*Eu;
    Uy = Rblk*Dy3*Eu;
    %
    % Director rotation beta = theta x n and its global derivatives.
    BetaG  = BetaOp*Er;
    BetaGx = Dx3*BetaG;
    BetaGy = Dy3*BetaG;
    Beta   = Rblk*BetaG;
    BetaX  = Rblk*BetaGx;
    BetaY  = Rblk*BetaGy;
    %
    r1 = 1:3:3*n_el;
    r2 = 2:3:3*n_el;
    r3 = 3:3:3*n_el;
    %
    % Membrane strains: symmetric tangential part of grad(u).
    E11 = Ux(r1,:);
    E22 = Uy(r2,:);
    G12 = Uy(r1,:) + Ux(r2,:);
    %
    % Changes of curvature from the tangential derivatives of beta.
    K11 = BetaX(r1,:);
    K22 = BetaY(r2,:);
    K12 = BetaY(r1,:) + BetaX(r2,:);
    %
    % Reissner-Mindlin transverse shear strains.
    S1 = Ux(r3,:) + Beta(r1,:);
    S2 = Uy(r3,:) + Beta(r2,:);
    %
    % Direct stiffness assembly in the global six-DOF coordinates.
    k_loc = zeros(6*n_el);
    k_loc = k_loc ...
        + sem2D.lame*C2*(E11.'*VD*E11 + E22.'*VD*E22) ...
        + sem2D.nu*sem2D.lame*C2*(E11.'*VD*E22 + E22.'*VD*E11) ...
        + sem2D.G*C2*(G12.'*VD*G12);
    %
    k_loc = k_loc ...
        + sem2D.lame*C1*(K11.'*VD*K11 + K22.'*VD*K22) ...
        + sem2D.nu*sem2D.lame*C1*(K11.'*VD*K22 + K22.'*VD*K11) ...
        + sem2D.G*C1*(K12.'*VD*K12);
    %
    k_loc = k_loc ...
        + k_sc*sem2D.G*C2*(S1.'*VD*S1 + S2.'*VD*S2);
    %
    % Consistent translational and rotary inertia in global coordinates.
    m_loc = zeros(6*n_el);
    M3 = kron(VD,eye(3));
    m_loc = m_loc ...
        + sem2D.rho*C2*(Eu.'*M3*Eu) ...
        + sem2D.rho*C1*(Beta(r1,:).'*VD*Beta(r1,:) ...
                       + Beta(r2,:).'*VD*Beta(r2,:));
    %
    % Drilling stabilization only; no artificial drilling inertia.
    th3 = NormalOp*Er;
    if isfield(sem2D,'etaD')
        etaD = sem2D.etaD;
    else
        etaD = 1e-6;
    end
    kd = etaD*sem2D.G*C2;
    k_loc = k_loc + kd*(th3.'*VD*th3);
    %
    % Numerical symmetry cleanup.
    k_loc = 0.5*(k_loc + k_loc.');
    m_loc = 0.5*(m_loc + m_loc.');
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
    VD = sem2D.VD * diag(sem2D.J(1,:,el));
    %
    QDxi_dxidx    = reshape(sem2D.InvJmat(1,1,:,el),n_el,1).*sem2D.Q1xi;
    QDxi_dxidy    = reshape(sem2D.InvJmat(2,1,:,el),n_el,1).*sem2D.Q1xi;
    QDeta_detadx  = reshape(sem2D.InvJmat(1,2,:,el),n_el,1).*sem2D.Q1eta;
    QDeta_detady  = reshape(sem2D.InvJmat(2,2,:,el),n_el,1).*sem2D.Q1eta;
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

    % k_loc = 0.5.*(k_loc + transpose(k_loc));
    % m_loc = 0.5.*(m_loc + transpose(m_loc));
end