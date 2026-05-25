clc; clear; close all;

addpath('../sem_core')
addpath('../sem_core/geometry')

% ============================================================
% INPUT GEOMETRY
% ============================================================
FileName = 'vertCyl_';
numPatch = 1;

Nurbs2D = iga2Dmesh(FileName,numPatch,1);

% ============================================================
% PARAMETRIC POINT
% ============================================================
u = 0;
v = 0;

% ============================================================
% DEGREE / KNOTS
% ============================================================
p = Nurbs2D.order{1}(1) - 1;
q = Nurbs2D.order{1}(2) - 1;

U = Nurbs2D.knots.U{1};
V = Nurbs2D.knots.V{1};

nu = Nurbs2D.number{1}(1);
nv = Nurbs2D.number{1}(2);

CP = Nurbs2D.cPoints{1};

% ============================================================
% FIND SPAN
% ============================================================
iu = findspan(nu, p, u, U);
iv = findspan(nv, q, v, V);

Nu = dersbasisfuns(iu, u, p, 2, U);
Nv = dersbasisfuns(iv, v, q, 2, V);

CPw = CP(:, iu-p:iu, iv-q:iv);

W  = squeeze(CPw(4,:,:));
Xw = squeeze(CPw(1,:,:)) .* W;
Yw = squeeze(CPw(2,:,:)) .* W;
Zw = squeeze(CPw(3,:,:)) .* W;

% ============================================================
% BASIS TENSORS
% ============================================================
Nu0 = Nu(1,:); Nu1 = Nu(2,:); Nu2 = Nu(3,:);
Nv0 = Nv(1,:); Nv1 = Nv(2,:); Nv2 = Nv(3,:);

N    = Nu0' * Nv0;
NuN  = Nu1' * Nv0;
NvN  = Nu0' * Nv1;

NuuN = Nu2' * Nv0;
NuvN = Nu1' * Nv1;
NvvN = Nu0' * Nv2;

% ============================================================
% WEIGHTED SUMS
% ============================================================
w   = sum(sum(N    .* W));
wu  = sum(sum(NuN  .* W));
wv  = sum(sum(NvN  .* W));

wuu = sum(sum(NuuN .* W));
wuv = sum(sum(NuvN .* W));
wvv = sum(sum(NvvN .* W));

X  = sum(sum(N    .* Xw));
Y  = sum(sum(N    .* Yw));
Z  = sum(sum(N    .* Zw));

Xu = sum(sum(NuN  .* Xw));
Yu = sum(sum(NuN  .* Yw));
Zu = sum(sum(NuN  .* Zw));

Xv = sum(sum(NvN  .* Xw));
Yv = sum(sum(NvN  .* Yw));
Zv = sum(sum(NvN  .* Zw));

Xuu = sum(sum(NuuN .* Xw));
Yuu = sum(sum(NuuN .* Yw));
Zuu = sum(sum(NuuN .* Zw));

Xuv = sum(sum(NuvN .* Xw));
Yuv = sum(sum(NuvN .* Yw));
Zuv = sum(sum(NuvN .* Zw));

Xvv = sum(sum(NvvN .* Xw));
Yvv = sum(sum(NvvN .* Yw));
Zvv = sum(sum(NvvN .* Zw));

% ============================================================
% QUOTIENT RULE (STABLE FORM)
% ============================================================
S = [X; Y; Z] / w
norm(S)

Su = ([Xu;Yu;Zu]*w - [X;Y;Z]*wu) / w^2
norm(Su)
Sv = ([Xv;Yv;Zv]*w - [X;Y;Z]*wv) / w^2
norm(Sv)

Suu = ([Xuu;Yuu;Zuu]*w^2 ...
      - 2*[Xu;Yu;Zu]*w*wu ...
      - [X;Y;Z]*(w*wuu - 2*wu^2)) / w^3;

Suv = ([Xuv;Yuv;Zuv]*w^2 ...
      - [Xu;Yu;Zu]*w*wv ...
      - [Xv;Yv;Zv]*w*wu ...
      - [X;Y;Z]*(w*wuv - 2*wu*wv)) / w^3;

Svv = ([Xvv;Yvv;Zvv]*w^2 ...
      - 2*[Xv;Yv;Zv]*w*wv ...
      - [X;Y;Z]*(w*wvv - 2*wv^2)) / w^3;

% ============================================================
% FUNDAMENTAL FORMS
% ============================================================
n = cross(Su,Sv);
n = n / norm(n);

E = dot(Su,Su);
F = dot(Su,Sv);
G = dot(Sv,Sv);

L = dot(Suu,n);
M = dot(Suv,n);
Nf = dot(Svv,n);

den = E*G - F^2;

% ============================================================
% CURVATURE
% ============================================================
K = (L*Nf - M^2) / den;

H = (E*Nf - 2*F*M + G*L) / (2*den);

disc = sqrt(max(H^2 - K, 0));

k1 = H + disc;
k2 = H - disc;

% ============================================================
% OUTPUT
% ============================================================
fprintf('\n--- CURVATURE RESULTS ---\n');
fprintf('Point S  = [%f %f %f]\n', S(1),S(2),S(3));
fprintf('K        = %.10f\n', K);
fprintf('H        = %.10f\n', H);
fprintf('k1       = %.10f\n', k1);
fprintf('k2       = %.10f\n', k2);