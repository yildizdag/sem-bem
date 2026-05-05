function [INC,IEN,nel,nnp,nen] = connectivity(order,number)

% --------------------------------------------
%   Create connectivity arrays of a d=1,2 or 3 
%   dimensional patch
%   
%   INPUT:
%   order  - polynomial orders in each parametric direction
%   number - number of basis functions in each parametric direction
%
%   OUTPUT:
%   INC - NURBS co-ordinate array, (nnp x d)
%   IEN - element node array, (nen x nel)
%   nel - number of elements
%   nnp - number of global basis
%   nen - number of local basis
% --------------------------------------------

e = 0; 
A = 0;
if numel(order) == 1
    
    d = 1;
    p = order-1;
    n = number;
    nel = n-p;
    nnp = n;
    nen = p+1;
    INC = zeros(nnp,d);     
    IEN = zeros(nen,nel);   
    
    for i = 1:n
        A = A + 1;
        INC(A,1) = i;
        
        if i >= p+1
            e = e + 1;
            
            for iloc = 0:p
                B = A - iloc;
                b = iloc + 1;
                
                IEN(b,e) = B;
            end
        end
    end
                
elseif numel(order) == 2
    
    d = 2;
    p = order(1)-1; q = order(2)-1;
    n = number(1);  m = number(2);
    nel = (n-p)*(m-q);
    nnp = n*m;
    nen = (p+1)*(q+1);
    INC = zeros(nnp,d);     
    IEN = zeros(nen,nel);   
    
    for j = 1:m
        for i = 1:n
            A = A + 1;
            INC(A,1) = i;
            INC(A,2) = j;
        
            if i >= p+1 && j >= q+1
                e = e + 1;
                
                for jloc = 0:q
                    for iloc = 0:p
                        B = A - jloc*n - iloc;
                        b = jloc*(p+1) + iloc + 1;

                        IEN(b,e) = B;
                    end
                end 
            end      
        end
    end
    
elseif numel(order) == 3
    
    d = 3;
    p = order(1)-1; q = order(2)-1; r = order(3)-1;
    n = number(1);  m = number(2);  l = number(3);
    nel = (n-p)*(m-q)*(l-r);
    nnp = n*m*l;
    nen = (p+1)*(q+1)*(r+1);
    INC = zeros(nnp,d);     
    IEN = zeros(nen,nel);   
    
    for k = 1:l
        for j = 1:m
            for i = 1:n
                A = A + 1;
                INC(A,1) = i;
                INC(A,2) = j;
                INC(A,3) = k;

                if i >= p+1 && j >= q+1 && k >= r+1
                    e = e + 1;
                    
                    for kloc = 0:r
                        for jloc = 0:q
                            for iloc = 0:p
                                B = A - kloc*m*n - jloc*n - iloc;
                                b = kloc*(p+1)*(q+1) + jloc*(p+1) + iloc + 1;
                                IEN(b,e) = B;
                            end
                        end
                    end 

                end

            end
        end
    end
    
end