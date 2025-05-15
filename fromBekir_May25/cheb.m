function [FT,IT] = cheb (space)

% ------------ Forward and Backward Transformation Matrix -----------------
% This code calculates forward and backward transformation of Tchebychev 
% polynomials. An analytical function can be expressed using Tchebychev
% polynomials as,
%
%                   y(x) = a0*T0+a1*T1+a2*T3+...
%
% Above function can also be expressed as,
% 
%                   y = IT * a
% or
%
%                   a = FT * y
%
% where   y  : sampled function
%         a  : coefficients of chebychev polynomials
%         IT : backward (inverse) transformation
%         FT : forward transformation
%
%INPUTS:
%   - space   : space array of the interested domain
%
%OUTPUTS:
%   - IT      : backward (inverse) transformation
%   - FT      : forward transformation
%
% -----------------------------------------------------------------

IT = zeros(space.N);
FT = zeros(space.N);

for i=1:space.N
    for j=1:space.N
            IT(i,j) = cos((j-1)*pi*(space.N-i)/(space.N-1));
    end
end

FT = IT^-1;