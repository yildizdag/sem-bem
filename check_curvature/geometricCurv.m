clc; clear; close all;

addpath('../sem_core')
addpath('../sem_core/geometry')

% ============================================================
% LOAD SURFACE
% ============================================================
FileName = 'vertCyl_';
Nurbs2D = iga2Dmesh(FileName,1,1);

u = 0;
v = 0;

h = 1e-6;

% ============================================================
% SURFACE EVALUATION FUNCTION (local helper)
% ============================================================
evalSurf = @(u,v) evalPoint(Nurbs2D,u,v);

% ============================================================
% BASE POINT
% ============================================================
[S, Su, Sv] = evalSurf(u,v);

n0 = cross(Su,Sv);
n0 = n0 / norm(n0);

E = dot(Su,Su);
F = dot(Su,Sv);
G = dot(Sv,Sv);

I = [E F; F G];

% ============================================================
% NORMAL VARIATION (ONLY FIRST DERIVATIVES USED)
% ============================================================
[S1, Su1, Sv1] = evalSurf(u+h,v);
n_u = cross(Su1,Sv1); n_u = n_u / norm(n_u);

[S2, Su2, Sv2] = evalSurf(u,v+h);
n_v = cross(Su2,Sv2); n_v = n_v / norm(n_v);

dn_u = (n_u - n0)/h;
dn_v = (n_v - n0)/h;

% ============================================================
% PROJECT ON TANGENT BASIS
% ============================================================
B = [dot(dn_u,Su) dot(dn_u,Sv);
     dot(dn_v,Su) dot(dn_v,Sv)];

% ============================================================
% SHAPE OPERATOR
% ============================================================
Sop = I \ B;

kappa = eig(Sop);

k1 = max(kappa);
k2 = min(kappa);

K = det(Sop);
H = trace(Sop)/2;

% ============================================================
% OUTPUT
% ============================================================
fprintf('\n--- FIRST-DERIVATIVE CURVATURE ---\n');
fprintf('K  = %.10f\n', K);
fprintf('H  = %.10f\n', H);
fprintf('k1 = %.10f\n', k1);
fprintf('k2 = %.10f\n', k2);

% ============================================================
% LOCAL FUNCTION
% ============================================================
function [S, Su, Sv] = evalPoint(Nurbs2D,u,v)

p = Nurbs2D.order{1}(1)-1;
q = Nurbs2D.order{1}(2)-1;

U = Nurbs2D.knots.U{1};
V = Nurbs2D.knots.V{1};

nu = Nurbs2D.number{1}(1);
nv = Nurbs2D.number{1}(2);

CP = Nurbs2D.cPoints{1};

iu = findspan(nu,p,u,U);
iv = findspan(nv,q,v,V);

Nu = dersbasisfuns(iu,u,p,1,U);
Nv = dersbasisfuns(iv,v,q,1,V);

CPw = CP(:,iu-p:iu,iv-q:iv);

W  = squeeze(CPw(4,:,:));
Xw = squeeze(CPw(1,:,:)).*W;
Yw = squeeze(CPw(2,:,:)).*W;
Zw = squeeze(CPw(3,:,:)).*W;

Nu0 = Nu(1,:); Nu1 = Nu(2,:);
Nv0 = Nv(1,:); Nv1 = Nv(2,:);

N  = Nu0'*Nv0;
NuN = Nu1'*Nv0;
NvN = Nu0'*Nv1;

w  = sum(sum(N.*W));
wu = sum(sum(NuN.*W));
wv = sum(sum(NvN.*W));

X  = sum(sum(N.*Xw));
Y  = sum(sum(N.*Yw));
Z  = sum(sum(N.*Zw));

Xu = sum(sum(NuN.*Xw));
Yu = sum(sum(NuN.*Yw));
Zu = sum(sum(NuN.*Zw));

Xv = sum(sum(NvN.*Xw));
Yv = sum(sum(NvN.*Yw));
Zv = sum(sum(NvN.*Zw));

S  = [X;Y;Z]/w;
Su = ([Xu;Yu;Zu]*w - [X;Y;Z]*wu)/w^2;
Sv = ([Xv;Yv;Zv]*w - [X;Y;Z]*wv)/w^2;

end