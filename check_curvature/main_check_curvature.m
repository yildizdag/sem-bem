%==========================================================================
% Project       : TÜBİTAK 3501 (125M858)
%
% Description   :
%  This code checks curvature computation using NURBS surfaces.
%
% Authors       : M. Erden Yildizdag, Bekir Bediz
%==========================================================================
clc; clear; close all;
addpath('../sem_core')
addpath('../sem_core/geometry')
%-Read the Geometry:
FileName = 'vertCyl_';
numPatch = 1; %Enter #Patches
%
Nurbs2D = iga2Dmesh(FileName,numPatch,1);
%
% Nurbs2D = hrefine2D(Nurbs2D,1,0,0);
%
figure;
iga2DmeshPlotNURBS(Nurbs2D);
%
np_u = 5;

np_v = 5;
epsilon = 1E-6;
%
u_sample = 2;
v_sample = 3;
%
k=1;
%
iu = findspan(Nurbs2D.number{1}(1),Nurbs2D.order{1}(1)-1,u_sample,Nurbs2D.knots.U{1});
iv = findspan(Nurbs2D.number{1}(2),Nurbs2D.order{1}(2)-1,v_sample,Nurbs2D.knots.V{1});
%
dNu = dersbasisfuns(iu,u_sample,Nurbs2D.order{1}(1)-1,2,Nurbs2D.knots.U{1});
dNv = dersbasisfuns(iv,v_sample,Nurbs2D.order{1}(2)-1,2,Nurbs2D.knots.V{1});
CP = Nurbs2D.cPoints{1}(:,iu-Nurbs2D.order{1}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
[dR,dS] = derRat2DBasisFuns(dNu,dNv,Nurbs2D.order{1}(1),Nurbs2D.order{k}(2),CP,2,2);
S0 = (dS(:,1,1))
r = sqrt(S0(1)^2+S0(2)^2)
A1 = dS(:,2,1); A2 = dS(:,1,2);
SN = cross(A1,A2)/norm(cross(A1,A2))
% F1 = [dot(A1,A1), dot(A1,A2)
%     dot(A1,A2), dot(A2,A2)];
% F2 = [dot(dS(:,3,1),N), dot(dS(:,2,2),N)
%     dot(dS(:,2,2),N), dot(dS(:,1,3),N)];
% F = F1\F2;
%
% First fundamental form
E = dot(A1,A1);
F = dot(A1,A2);
G = dot(A2,A2);
% Second fundamental form
L = dot(dS(:,3,1),SN);
M = dot(dS(:,2,2),SN);
N = dot(dS(:,1,3),SN);
%
den = E*G - F^2;
%
K = (L*N - M^2)/den;
%
H = (E*N - 2*F*M + G*L)/(2*den);
%
% kappa = 1E-4.*round(eig(F)./1E-4)
k1 = H + sqrt(H^2 - K)
k2 = H - sqrt(H^2 - K)
norm(dS(:,3,1))
