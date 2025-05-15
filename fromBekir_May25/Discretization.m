function [FT,IT,D,s,V,Q1,Q2,space] = Discretization (L, polynum, axis)

% ------------ DISCRETIZATION OF COORDINATES ----------------------
% This code discritize the geometry for the given direction and 
% calculates the necessary parameters for Tchebychev solution:
%   - Forward and backward transformation [FT,IT]
%   - Derivative matrix [D]
%   - Definite integral vector [vint]
%   - Discretization points [s]
%   - Inner product matrix [V]
%   - Spatial derivative matrix [Q]
%   - Defining the space in each axis [space]
%
%INPUTS:
%   - L       : discreatization length
%   - polynum : number of polynomials used in discretization 
%               (# of discretization points) 
%   - axis    : axis in which the discretization will be performed (since 
%               the axis (z axis) is different than the other two axes)
%
%OUTPUTS:
%   - FT    : forward transformation matrix [a=FT*y]
%   - IT    : inverse transformation matrix [y=IT*a]
%   - D     : spatial differentiation matrix [b=Da]
%   - s     : physical sampling points in the domain
%   - V     : inner product matrix [int(f*g)=f_t*V*g]
%   - Q1    : spatial derivative matrix [Q=IT*D*FT]
%   - space : space array for each domain
%
% -----------------------------------------------------------------

spn=polynum;    % number of Tcheb polynomials

% defining the space (coordinates)
switch axis
    case {('xi'),('x')}
        space.a=-L/2;   
        space.b=L/2;  
        space.N=spn;
    case {('eta'), ('y')}
        space.a=-L/2;   
        space.b=L/2;  
        space.N=spn;
    case {('zeta'), ('z')}
        space.a=0;   
        space.b=L;  
        space.N=spn;
end

[FT,IT]=cheb(space);
[D] = derivative(space);
s=slobat(space); 
V=InnerProduct(space);
Q1 = IT*D*FT;
Q2 = IT*(D^2)*FT;

end

