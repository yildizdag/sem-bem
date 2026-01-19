function [FT,IT] = transform(N)
IT = zeros(N);
for i=1:N
    for j=1:N
        IT(i,j) = cos((j-1)*pi*(N-i)/(N-1));
    end
end
FT = IT^-1;