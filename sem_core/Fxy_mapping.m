function [Fxy] = Fxy_mapping(spnx, spny, Fx, Fy)

% ---------------------------------------------------------------
% This code calculates the Fxyz (3D vector mapping) from Fx, Fy, and Fz.
%
%INPUTS:
%   - spnx       : number of sampling points in x direction
%   - spny       : number of sampling points in y direction
%   - spnz       : number of sampling points in z direction
%   - Fx         : function describing the x direction
%   - Fy         : function describing the y direction
%   - Fz         : function describing the z direction
%s
%OUTPUTS:
%   - Fxyz       : function describing the whole xyz domain
%  
% ---------------------------------------------------------------

spn_xy = spnx*spny;

Fxy = zeros(spn_xy,spn_xy);

for j=1:spnx
    for k=1:spny        
        for j2=1:spnx
            for k2=1:spny
                count = (j-1)*spny+k;
                count2 = (j2-1)*spny+k2;
                Fxy(count2,count) = Fxy(count2,count) + Fx(j2,j)*Fy(k2,k);
            end
        end
    end
end
