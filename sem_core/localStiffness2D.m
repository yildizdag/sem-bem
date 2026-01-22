function k_loc = localStiffness2D(sem2D,el)
%
ldof = sem2D.local_dof;
nconn = sem2D.conn(el,ldof:ldof:end)./ldof;

if sem2D.ET == 1
    %
%     BT = kron(sem2D.BT_xi,sem2D.BT_eta);
%     FT = BT^-1;
%     Dxi = kron(sem2D.D_xi,eye(sem2D.N,sem2D.N));
%     Deta = kron(sem2D.D_eta,eye(sem2D.N,sem2D.N));
    Q1xi = kron(eye(sem2D.N,sem2D.N),sem2D.Q1_xi);
    Q1eta = kron(eye(sem2D.N,sem2D.N),sem2D.Q1_eta);
    dx_dxi = Q1xi*sem2D.nodes(nconn,1);
    dy_dxi = Q1xi*sem2D.nodes(nconn,2);
    dx_deta = Q1eta*sem2D.nodes(nconn,1);
    dy_deta = Q1eta*sem2D.nodes(nconn,2);
    %
    J = dy_deta.*dx_dxi-dy_dxi.*dx_deta
    %J = (sem2D.Q1)*sem2D.nodes(nconn,1);
    k_loc = zeros(sem2D.N*sem2D.N);
elseif sem2D.ET == 2
    k_loc = zeros(sem1D.local_dof*sem1D.N);
    J = (sem1D.Q1)*sem1D.nodes(nconn,1);
    k_loc(2:2:end,2:2:end) = k_loc(2:2:end,2:2:end) + (sem1D.EI).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(1:2:end,1:2:end) = k_loc(1:2:end,1:2:end) + (sem1D.kGA).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(1:2:end,2:2:end) = k_loc(1:2:end,2:2:end) - (sem1D.kGA).*((sem1D.Q1./J)'*sem1D.V).*J;
    k_loc(2:2:end,1:2:end) = k_loc(2:2:end,1:2:end) - (sem1D.kGA).*(sem1D.V*(sem1D.Q1./J)).*J;
    k_loc(2:2:end,2:2:end) = k_loc(2:2:end,2:2:end) + (sem1D.kGA).*sem1D.V.*J;
end