function [D] = derivative(space)

% ---------------- Spatial Derivative Matrix ------------------------------
% This code calclates the spatial derivative matrix of a function expressed
% using Tchebychev polynomials. An analytical function  and its spatial 
% derivative can be expressed using Tchebychev polynomials as,
%
%                   y(x) = a0*T0+a1*T1+a2*T3+...
%
%                   y'(x) = b0*T0+b1*T1+b2*T3+...   
%
% Therefore, a spatial derivative matrix can be defined between the
% coefficients, a and b as,
%
%                  b = D*a
%
% Note that the spatial differentiation matrix is calculated for the domain
% (-1,1). Thus, in order to transform it into the domain (l1,l2), a type of
% coordinate transformation is needed. Let,
%
%                 e(x) = 2/(l2-l1)*x-(l2+l1)/(l2-l1),     e ---> (-1,1)
%                                                         x ---> (l1,l2)
%
% Taking the derivative of the above equation,
%
%                 de/dx = 2/(l2-l1) * x
%
% Therefore, the spatial derivative of the domain of interest can be as,
%
%                 D(for (l1,l2)) = D (for (-1,1)) / scale
%
% where  scale = (l2-l1)/2
%        D (for (-1,1)) = de/dx
%        D(for (l1,l2)) = dx
%
%INPUTS:
%   - space   : space array of the interested domain
%
%OUTPUTS:
%   - D       : Spatial derivative matrix 
%
% -------------------------------------------------------------------------



%defining the scale factor
scale =(space.b-space.a)/2;

if rem(space.N,2) == 0
    evenorodd=1;
else
    evenorodd=0;
end

D = zeros(space.N,space.N);

if evenorodd == 1
    DN = (space.N)/2;
    for i=0:DN-1
        D(1,2*i+2) = 2*i+1;
    end
    
    for i=2:space.N
        if rem(i,2) == 0
            for j=1:DN-1
                D(i,2*j+1) = 4*j;
            end
        else
            for j=1:DN-1
                D(i,2*j+2) = 2*(2*j+1);
            end
        end
    end
    
else
    
    DN = (space.N-1)/2;
    for i=0:DN-1
        D(1,2*i+2) = 2*i+1;
    end
    
    for i=2:space.N
        if rem(i,2) == 0
            for j=1:DN
                D(i,2*j+1) = 4*j;
            end
        else
            for j=1:DN-1
                D(i,2*j+2) = 2*(2*j+1);
            end
        end
    end
    
    
end

D=triu(D);

D=D/scale;