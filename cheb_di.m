function [v_di] = cheb_di(space)
% This gives the vector vecint such that when you dot
% the vector of Tchebychev coefficents with 
% vecint you get the definite integral of the function on the
% interval (xleft,xright).

v_di = zeros(space.N,1);

for i = 1:2:space.N
    v_di(i) = (space.b-space.a)/(1-(i-1)^2);
end