function [FT,IT] = cheb(space)

IT = zeros(space.N);

for i=1:space.N
    for j=1:space.N
        IT(i,j) = cos((j-1)*pi*(space.N-i)/(space.N-1));
    end
end

FT = IT^-1;