function [K_2D, M_2D] = Mass_and_Stiffness_FSDT_DC(polynum_xi, polynum_eta, ...
                   ro, lame, mu, pois, k_sc, h,...
                   VD, QDx, QDy, beta_1, beta_2)

% ------------------------------------------------------

spnx = polynum_xi;
spny = polynum_eta;
spnxy = spnx*spny; 

ub = 1;
ue = spnxy;
vb = spnxy+1;
ve = 2*spnxy;
wb = 2*spnxy+1;
we = 3*spnxy;
phi_xb = 3*spnxy+1;
phi_xe = 4*spnxy;
phi_yb = 4*spnxy+1;
phi_ye = 5*spnxy;

% -------------------------------
M_2D = zeros(5*spnxy, 5*spnxy);
M_2D(ub:ue, ub:ue) = ro*h*VD;
M_2D(vb:ve, vb:ve) = ro*h*VD;
M_2D(wb:we, wb:we) = ro*h*VD;
M_2D(phi_xb:phi_xe, phi_xb:phi_xe) = ro*h^3/12*VD;    
M_2D(phi_yb:phi_ye, phi_yb:phi_ye) = ro*h^3/12*VD; 
% -------------------------------


% -------------------------------
C1 = 1/12*h^3;
C2 = h;
K_2D = zeros(5*spnxy, 5*spnxy);     
% eps_x*sigma_x
K_2D(ub:ue,ub:ue) = K_2D(ub:ue,ub:ue) + lame*C2*QDx'*VD*QDx;
K_2D(ub:ue,wb:we) = K_2D(ub:ue,wb:we) + lame*C2*beta_1*QDx'*VD;
K_2D(wb:we,ub:ue) = K_2D(wb:we,ub:ue) + lame*C2*beta_1*VD*QDx;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + lame*C2*beta_1^2*VD;
K_2D(phi_xb:phi_xe,phi_xb:phi_xe) = K_2D(phi_xb:phi_xe,phi_xb:phi_xe) + lame*C1*QDx'*VD*QDx;

K_2D(ub:ue,vb:ve) = K_2D(ub:ue,vb:ve) + pois*lame*C2*QDx'*VD*QDy;
K_2D(ub:ue,wb:we) = K_2D(ub:ue,wb:we) + pois*lame*C2*beta_2*QDx'*VD;
K_2D(wb:we,vb:ve) = K_2D(wb:we,vb:ve) + pois*lame*C2*beta_1*VD*QDy;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + pois*lame*C2*beta_1*beta_2*VD;
K_2D(phi_xb:phi_xe,phi_yb:phi_ye) = K_2D(phi_xb:phi_xe,phi_yb:phi_ye) + pois*lame*C1*QDx'*VD*QDy;

% eps_y*sigma_y
K_2D(vb:ve,ub:ue) = K_2D(vb:ve,ub:ue) + pois*lame*C2*QDy'*VD*QDx;
K_2D(vb:ve,wb:we) = K_2D(vb:ve,wb:we) + pois*lame*C2*beta_1*QDy'*VD;
K_2D(wb:we,ub:ue) = K_2D(wb:we,ub:ue) + pois*lame*C2*beta_2*VD*QDx;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + pois*lame*C2*beta_1*beta_2*VD;
K_2D(phi_yb:phi_ye,phi_xb:phi_xe) = K_2D(phi_yb:phi_ye,phi_xb:phi_xe) + pois*lame*C1*QDy'*VD*QDx;

K_2D(vb:ve,vb:ve) = K_2D(vb:ve,vb:ve) + lame*C2*QDy'*VD*QDy;
K_2D(vb:ve,wb:we) = K_2D(vb:ve,wb:we) + lame*C2*beta_2*QDy'*VD;
K_2D(wb:we,vb:ve) = K_2D(wb:we,vb:ve) + lame*C2*beta_2*VD*QDy;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + lame*C2*beta_2^2*VD;
K_2D(phi_yb:phi_ye,phi_yb:phi_ye) = K_2D(phi_yb:phi_ye,phi_yb:phi_ye) + lame*C1*QDy'*VD*QDy;

% gamma_xy*tau_xy
K_2D(ub:ue,ub:ue) = K_2D(ub:ue,ub:ue) + mu*C2*QDy'*VD*QDy;
K_2D(ub:ue,vb:ve) = K_2D(ub:ue,vb:ve) + mu*C2*QDy'*VD*QDx;
K_2D(vb:ve,ub:ue) = K_2D(vb:ve,ub:ue) + mu*C2*QDx'*VD*QDy;
K_2D(vb:ve,vb:ve) = K_2D(vb:ve,vb:ve) + mu*C2*QDx'*VD*QDx;
K_2D(phi_xb:phi_xe,phi_xb:phi_xe) = K_2D(phi_xb:phi_xe,phi_xb:phi_xe) + mu*C1*QDy'*VD*QDy;
K_2D(phi_xb:phi_xe,phi_yb:phi_ye) = K_2D(phi_xb:phi_xe,phi_yb:phi_ye) + mu*C1*QDy'*VD*QDx;
K_2D(phi_yb:phi_ye,phi_xb:phi_xe) = K_2D(phi_yb:phi_ye,phi_xb:phi_xe) + mu*C1*QDx'*VD*QDy;
K_2D(phi_yb:phi_ye,phi_yb:phi_ye) = K_2D(phi_yb:phi_ye,phi_yb:phi_ye) + mu*C1*QDx'*VD*QDx;

% gamma_xz*tau_xz
K_2D(ub:ue,ub:ue) = K_2D(ub:ue,ub:ue) + k_sc^2*mu*C2*beta_1^2*VD;
K_2D(ub:ue,wb:we) = K_2D(ub:ue,wb:we) - k_sc^2*mu*C2*beta_1*VD*QDx;
K_2D(ub:ue,phi_xb:phi_xe) = K_2D(ub:ue,phi_xb:phi_xe) - k_sc^2*mu*C2*beta_1*VD;
K_2D(wb:we,ub:ue) = K_2D(wb:we,ub:ue) - k_sc^2*mu*C2*beta_1*QDx'*VD;
K_2D(phi_xb:phi_xe,ub:ue) = K_2D(phi_xb:phi_xe,ub:ue) - k_sc^2*mu*C2*beta_1*VD;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + k_sc^2*mu*C2*QDx'*VD*QDx;
K_2D(wb:we,phi_xb:phi_xe) = K_2D(wb:we,phi_xb:phi_xe) + k_sc^2*mu*C2*QDx'*VD;
K_2D(phi_xb:phi_xe,wb:we) = K_2D(phi_xb:phi_xe,wb:we) + k_sc^2*mu*C2*VD*QDx;
K_2D(phi_xb:phi_xe,phi_xb:phi_xe) = K_2D(phi_xb:phi_xe,phi_xb:phi_xe) + k_sc^2*mu*C2*VD;

% gamma_yz*tau_yz
K_2D(vb:ve,vb:ve) = K_2D(vb:ve,vb:ve) + k_sc^2*mu*C2*beta_2*beta_2*VD;
K_2D(vb:ve,wb:we) = K_2D(vb:ve,wb:we) - k_sc^2*mu*C2*beta_2*VD*QDy;
K_2D(vb:ve,phi_yb:phi_ye) = K_2D(vb:ve,phi_yb:phi_ye) - k_sc^2*mu*C2*beta_2*VD;
K_2D(wb:we,vb:ve) = K_2D(wb:we,vb:ve) - k_sc^2*mu*C2*beta_2*QDy'*VD;
K_2D(phi_yb:phi_ye,vb:ve) = K_2D(phi_yb:phi_ye,vb:ve) - k_sc^2*mu*C2*beta_2*VD;
K_2D(wb:we,wb:we) = K_2D(wb:we,wb:we) + k_sc^2*mu*C2*QDy'*VD*QDy;
K_2D(wb:we,phi_yb:phi_ye) = K_2D(wb:we,phi_yb:phi_ye) + k_sc^2*mu*C2*QDy'*VD;
K_2D(phi_yb:phi_ye,wb:we) = K_2D(phi_yb:phi_ye,wb:we) + k_sc^2*mu*C2*VD*QDy;
K_2D(phi_yb:phi_ye,phi_yb:phi_ye) = K_2D(phi_yb:phi_ye,phi_yb:phi_ye) + k_sc^2*mu*C2*VD;
% -------------------------------



end
