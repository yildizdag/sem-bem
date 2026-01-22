function k_loc = localStiffness2D(sem2D,el)
%
ldof = sem2D.local_dof;
nconn = sem2D.conn(el,ldof:ldof:end)./ldof;
%
if sem2D.ET == 1
    %
    VD = sem2D.VD * diag(sem2D.J(nconn));
    %
    QDxi_dxidx    = diag(reshape(sem2D.InvJmat(1,1,nconn),sem2D.N*sem2D.N,1))*sem2D.Q1xi;
    QDxi_dxidy    = diag(reshape(sem2D.InvJmat(1,2,nconn),sem2D.N*sem2D.N,1))*sem2D.Q1xi;
    QDeta_detadx  = diag(reshape(sem2D.InvJmat(2,1,nconn),sem2D.N*sem2D.N,1))*sem2D.Q1eta;
    QDeta_detady  = diag(reshape(sem2D.InvJmat(1,2,nconn),sem2D.N*sem2D.N,1))*sem2D.Q1eta;
    %
    QDx = QDxi_dxidx + QDeta_detadx;
    QDy = QDxi_dxidy + QDeta_detady;
    %
    k_loc = zeros(sem2D.N*sem2D.N);
end