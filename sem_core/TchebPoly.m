function [Tcheb] = TchebPoly(degree,fine)

Tcheb = zeros(length(fine),degree+1);

Tcheb(:,1) = ones(length(fine),1);
Tcheb(:,2) = fine;

for i = 3:degree+1
        Tcheb(:,i) = 2*Tcheb(:,2).*Tcheb(:,i-1)-Tcheb(:,i-2);    
end
