function [dR,dC] = derRat1DBasisFuns(dNu,o,CP,du)
    w = CP(4,:)';
    wders = dNu*w;
    dR = zeros(o,du+1);
    for k = 0:du
        v = dNu(k+1,:)'.*w;
        for i = 0:k-1
            v = v - nchoosek(k,i)*(wders(k-i+1)*dR(:,i+1));
        end
        dR(:,k+1) = v/wders(1);
    end
    dC = zeros(3,du+1);
    for i = 1:du+1
            dC(:,i) = CP(1:3,:)*dR(:,i);
    end
end