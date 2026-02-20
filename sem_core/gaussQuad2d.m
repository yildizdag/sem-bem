function [xgp,wgp,ngp] = gaussQuad2d(ngpx,ngpy)
%--------------------------------------------------------------------------
%   Returns the points used for 2d Gaussian Quadrature
%   
%   INPUT:
%   ngpx - number of quadrature points along xi
%   ngpy - number of quadrature points along eta
%
%   OUTPUT:
%   xgp - Gauss points, (ngpx*ngpy x 2)
%   wgp - corresponding weights, (ngpx*ngpy x 1)
%   ngp - number of quadrature points (ngpx*ngpy)
%--------------------------------------------------------------------------
[xgp_xi,wgp_xi] = gaussQuad1d(ngpx);
[xgp_eta,wgp_eta] = gaussQuad1d(ngpy);
ngp = ngpx*ngpy;
xgp = zeros(ngp,2);
wgp = zeros(ngp,1);
for i = 1:ngpy
    k = (i-1)*ngpx;
    for j = 1:ngpx
        xgp(j+k,:) = [xgp_xi(j), xgp_eta(i)];
        wgp(j+k,1) = wgp_xi(j)*wgp_eta(i);
    end
end