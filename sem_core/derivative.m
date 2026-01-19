function D = derivative(space)

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

end