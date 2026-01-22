function [dR,dS] = derRat2DBasisFuns(dNu,dNv,o1,o2,CP,du,dv)

wij = reshape(CP(4,:,:),o1,o2);
dR = zeros(o1,o2,du+1,dv+1);
wders = zeros(du+1,dv+1);

for k = 0:du
    for l = 0:dv
        
        wders(k+1,l+1) = dNu(k+1,:) * wij * dNv(l+1,:)';
        
        v1 = dNu(k+1,:)' * dNv(l+1,:);
        v1 = v1 .* wij;
        
        for j = 1:l
            
            t1 = dR(:,:,k+1,l-j+1);
            v1 = v1 - (nchoosek(l,j)*wders(1,j+1)) .* t1;
            
        end
        
        for i = 1:k

            t1 = dR(:,:,k-i+1,l+1);
            v1 = v1 - (nchoosek(k,i) * wders(i+1,1)) .* t1;
            v2 = zeros((o1),(o2));
            
            for j=1:l
                
                t2 = dR(:,:,k-i+1,l-j+1);
                v2 = v2 + (nchoosek(l,j)*wders(i+1,j+1)) .* t2;
                
            end

            v1 = v1 - nchoosek(k,i).*v2;

        end
        
        dR(:,:,k+1,l+1) = v1 / wders(1,1);
         
    end
end

dS = zeros(3,du+1,dv+1);
for j = 1:du+1
    for k = 1:dv+1
        for i = 1:3
            dS(i,j,k) = reshape(dR(:,:,j,k),1,o1*o2) * reshape(CP(i,:,:),o1*o2,1);
        end
    end
end