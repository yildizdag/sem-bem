function [k_loc,m_loc] = local2D(sem2D,el)
%
n_dof = sem2D.local_dof;
nconn = sem2D.conn(el,n_dof:n_dof:end)./n_dof;
n_el = sem2D.N*sem2D.N;
%
if sem2D.ET == 1
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
elseif sem2D.ET == 2
    %
    k_loc1 = zeros(5*n_el);
    m_loc1 = zeros(5*n_el);
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
    k_loc1(1:5:end,1:5:end) = k_loc1(1:5:end,1:5:end) + sem2D.lame*sem2D.t*QDx'*VD*QDx;
    k_loc1(1:5:end,3:5:end) = k_loc1(1:5:end,3:5:end) + sem2D.lame*sem2D.t*beta_1*QDx'*VD;
    k_loc1(3:5:end,1:5:end) = k_loc1(3:5:end,1:5:end) + sem2D.lame*sem2D.t*beta_1*VD*QDx;
    k_loc1(3:5:end,3:5:end) = k_loc1(3:5:end,3:5:end) + sem2D.lame*sem2D.t*beta_1^2*VD;
    k_loc1(4:5:end,4:5:end) = k_loc1(4:5:end,4:5:end) + semd2D.lame*(sem2D.t^3/12)*QDx'*VD*QDx;
    %
    k_loc1(ub:ue,vb:ve) = k_loc1(ub:ue,vb:ve) + pois*lame*C2*QDx'*VD*QDy;
    k_loc1(ub:ue,wb:we) = k_loc1(ub:ue,wb:we) + pois*lame*C2*beta_2*QDx'*VD;
    k_loc1(wb:we,vb:ve) = k_loc1(wb:we,vb:ve) + pois*lame*C2*beta_1*VD*QDy;
    k_loc1(wb:we,wb:we) = k_loc1(wb:we,wb:we) + pois*lame*C2*beta_1*beta_2*VD;
    k_loc1(phi_xb:phi_xe,phi_yb:phi_ye) = k_loc1(phi_xb:phi_xe,phi_yb:phi_ye) + pois*lame*C1*QDx'*VD*QDy;

    % eps_y*sigma_y
    k_loc1(vb:ve,ub:ue) = k_loc1(vb:ve,ub:ue) + pois*lame*C2*QDy'*VD*QDx;
    k_loc1(vb:ve,wb:we) = k_loc1(vb:ve,wb:we) + pois*lame*C2*beta_1*QDy'*VD;
    k_loc1(wb:we,ub:ue) = k_loc1(wb:we,ub:ue) + pois*lame*C2*beta_2*VD*QDx;
    k_loc1(wb:we,wb:we) = k_loc1(wb:we,wb:we) + pois*lame*C2*beta_1*beta_2*VD;
    k_loc1(phi_yb:phi_ye,phi_xb:phi_xe) = k_loc1(phi_yb:phi_ye,phi_xb:phi_xe) + pois*lame*C1*QDy'*VD*QDx;

    k_loc1(vb:ve,vb:ve) = k_loc1(vb:ve,vb:ve) + lame*C2*QDy'*VD*QDy;
    k_loc1(vb:ve,wb:we) = k_loc1(vb:ve,wb:we) + lame*C2*beta_2*QDy'*VD;
    k_loc1(wb:we,vb:ve) = k_loc1(wb:we,vb:ve) + lame*C2*beta_2*VD*QDy;
    k_loc1(wb:we,wb:we) = k_loc1(wb:we,wb:we) + lame*C2*beta_2^2*VD;
    k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) = k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) + lame*C1*QDy'*VD*QDy;

    % gamma_xy*tau_xy
    k_loc1(ub:ue,ub:ue) = k_loc1(ub:ue,ub:ue) + mu*C2*QDy'*VD*QDy;
    k_loc1(ub:ue,vb:ve) = k_loc1(ub:ue,vb:ve) + mu*C2*QDy'*VD*QDx;
    k_loc1(vb:ve,ub:ue) = k_loc1(vb:ve,ub:ue) + mu*C2*QDx'*VD*QDy;
    k_loc1(vb:ve,vb:ve) = k_loc1(vb:ve,vb:ve) + mu*C2*QDx'*VD*QDx;
    k_loc1(phi_xb:phi_xe,phi_xb:phi_xe) = k_loc1(phi_xb:phi_xe,phi_xb:phi_xe) + mu*C1*QDy'*VD*QDy;
    k_loc1(phi_xb:phi_xe,phi_yb:phi_ye) = k_loc1(phi_xb:phi_xe,phi_yb:phi_ye) + mu*C1*QDy'*VD*QDx;
    k_loc1(phi_yb:phi_ye,phi_xb:phi_xe) = k_loc1(phi_yb:phi_ye,phi_xb:phi_xe) + mu*C1*QDx'*VD*QDy;
    k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) = k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) + mu*C1*QDx'*VD*QDx;

    % gamma_xz*tau_xz
    k_loc1(ub:ue,ub:ue) = k_loc1(ub:ue,ub:ue) + k_sc^2*mu*C2*beta_1^2*VD;
    k_loc1(ub:ue,wb:we) = k_loc1(ub:ue,wb:we) - k_sc^2*mu*C2*beta_1*VD*QDx;
    k_loc1(ub:ue,phi_xb:phi_xe) = k_loc1(ub:ue,phi_xb:phi_xe) - k_sc^2*mu*C2*beta_1*VD;
    k_loc1(wb:we,ub:ue) = k_loc1(wb:we,ub:ue) - k_sc^2*mu*C2*beta_1*QDx'*VD;
    k_loc1(phi_xb:phi_xe,ub:ue) = k_loc1(phi_xb:phi_xe,ub:ue) - k_sc^2*mu*C2*beta_1*VD;
    k_loc1(wb:we,wb:we) = k_loc1(wb:we,wb:we) + k_sc^2*mu*C2*QDx'*VD*QDx;
    k_loc1(wb:we,phi_xb:phi_xe) = k_loc1(wb:we,phi_xb:phi_xe) + k_sc^2*mu*C2*QDx'*VD;
    k_loc1(phi_xb:phi_xe,wb:we) = k_loc1(phi_xb:phi_xe,wb:we) + k_sc^2*mu*C2*VD*QDx;
    k_loc1(phi_xb:phi_xe,phi_xb:phi_xe) = k_loc1(phi_xb:phi_xe,phi_xb:phi_xe) + k_sc^2*mu*C2*VD;

    % gamma_yz*tau_yz
    k_loc1(vb:ve,vb:ve) = k_loc1(vb:ve,vb:ve) + k_sc^2*mu*C2*beta_2*beta_2*VD;
    k_loc1(vb:ve,wb:we) = k_loc1(vb:ve,wb:we) - k_sc^2*mu*C2*beta_2*VD*QDy;
    k_loc1(vb:ve,phi_yb:phi_ye) = k_loc1(vb:ve,phi_yb:phi_ye) - k_sc^2*mu*C2*beta_2*VD;
    k_loc1(wb:we,vb:ve) = k_loc1(wb:we,vb:ve) - k_sc^2*mu*C2*beta_2*QDy'*VD;
    k_loc1(phi_yb:phi_ye,vb:ve) = k_loc1(phi_yb:phi_ye,vb:ve) - k_sc^2*mu*C2*beta_2*VD;
    k_loc1(wb:we,wb:we) = k_loc1(wb:we,wb:we) + k_sc^2*mu*C2*QDy'*VD*QDy;
    k_loc1(wb:we,phi_yb:phi_ye) = k_loc1(wb:we,phi_yb:phi_ye) + k_sc^2*mu*C2*QDy'*VD;
    k_loc1(phi_yb:phi_ye,wb:we) = k_loc1(phi_yb:phi_ye,wb:we) + k_sc^2*mu*C2*VD*QDy;
    k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) = k_loc1(phi_yb:phi_ye,phi_yb:phi_ye) + k_sc^2*mu*C2*VD;
end