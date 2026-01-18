function [v_di] = cheb_di(space)

v_di = zeros(space.N,1);

for i = 1:2:space.N
    v_di(i) = (space.b-space.a)/(1-(i-1)^2);
end