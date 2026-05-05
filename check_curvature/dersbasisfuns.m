function ders = dersbasisfuns(i,u,p,k,U)

%--------------------------------------------------------------------------
%	Computes the non-zero basis functions, their derivatives and stores 
%   them in the array ders, where ders(k,j) denotes the kth derivative of
%   the function Ni-p+j,p. k=0,...,p and j=0,...,p.
%   Note: N(0)i,p gives the value of the basis function
%	Algorithm A2.3 from The NURBS book p.72
% 
%	INPUT: 
%   i - knot span index
%	u - parametric point
%	p - degree of the basis function
%	k - highest-degree derivative requested
%	U - knot vector
% 
%	OUTPUT:
%   ders - derivatives of the basis function, (p+1 x p+1) row-wise storage
%--------------------------------------------------------------------------

left  = zeros(1, p + 1);
right = zeros(1, p + 1);
ndu   = zeros(p + 1, p + 1);
a     = zeros(2, p + 1);
ders  = zeros(k + 1, p + 1);

ndu(1, 1) = 1;
for j = 1 : p
    left(j + 1)  = u - U(i + 1 - j);
    right(j + 1) = U(i + j) - u;
    saved = 0;
    for r = 0 : j - 1
        % lower triangle
        ndu(j + 1 , r + 1) = right(r + 2) + left(j - r + 1);
        temp = ndu(r + 1, j) / ndu(j + 1, r + 1);
        % upper triangle
        ndu(r + 1, j + 1) = saved + right(r + 2) * temp;
        saved = left(j - r + 1) * temp;
    end
    ndu(j + 1, j + 1) = saved; % ok
end

% load basis functions
for j = 1 : p+1    
    ders(1, j) = ndu(j, p + 1); % ok
end

% Compute the derivatives based on equation (2.10)
% loop over function index
for r = 0 : p
    % alternate rows in array a
    s1 = 1; % original 0
    s2 = 2; % original 1
    

    a(1, 1) = 1;
    % loop to compute hth derivative
    for h = 1 : k
        d = 0;
        rh = r - h;
        ph = p - h;
        
        if r >= h
            a(s2, 1) = a(s1, 1) / ndu(ph + 2, rh + 1);%
            d = a(s2, 1) * ndu(rh + 1, ph + 1);%
        end
        
        if rh >= -1
            j1 = 1;
        else
            j1 = -rh;
        end

        
        if r - 1 <= ph
            j2 = h - 1;
        else
            j2 = p - r;
        end
        
        for j = j1+1 : j2+1
            a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(ph + 2, rh + j);
            d = d + a(s2, j) * ndu(rh + j, ph + 1);
        end
        
        if r <= ph
            a(s2, h + 1) = -a(s1, h) / ndu(ph + 2, r + 1);
            d = d + a(s2, h + 1) * ndu(r + 1, ph + 1);           
        end
        
        ders(h + 1, r + 1) = d;
        % switch rows
        j  = s1;
        s1 = s2;
        s2 = j;
    end  
end

% multiply through by the correct factors of equation (2.10)
r = p;
for h = 1 : k
    for j = 0 : p
        ders(h + 1, j + 1) = ders(h + 1, j + 1) * r;
    end
    r = r * (p - h);
end