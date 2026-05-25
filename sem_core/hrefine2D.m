function Nurbs2D = hrefine2D(Nurbs2D,kk,ur,vr)
%----------------------------------------------------------------
% Adapted from The NURBS Book pp. 164-168
%----------------------------------------------------------------
if ur == 0
    disp(' ');
elseif ur > 0
    X = linspace(0,1,ur+2);
    X = X(2:end-1);
    UP = Nurbs2D.knots.U{kk};
    Pw = Nurbs2D.cPoints{kk};
    mp = numel(Nurbs2D.knots.U{kk}) - 1;
    p  = Nurbs2D.order{kk}(1)-1;
    np = mp - p - 1;
    r  = numel(X);
    %
    
    no1 = numel(UP) + numel(X);
    UQ  = zeros(1, no1);
    
    no2 = size(Pw,2) + numel(X);
    no3 = size(Pw,3);
    Qw  = zeros(4, no2, no3);
    
    a = findspan(np,p,X(1),UP) - 1; %* -1
    b = findspan(np,p,X(end),UP);   %* -1
    
    for i = 0 : a - p
        Qw(:,i + 1,:) = Pw(:,i + 1,:);
    end
    
    for i = b - 1 : np
        Qw(:,i + r + 1,:) = Pw(:,i + 1,:);
    end
    
    for i = 0 : a
        UQ(i + 1) = UP(i + 1);
    end
    
    for i = b + p : mp
        UQ(i + r + 1) = UP(i + 1);
    end
    
    i = b + p - 1 + 1;
    k = b + p + r;
    
    for j = r - 1 : -1 : 0
        while X(j + 1) <= UP(i) && i - 1 > a
            Qw(:,k - p - 1,:) = Pw(:,i - p - 1,:);
            UQ(k) = UP(i);
            k = k - 1;
            i = i - 1;
        end
        
        Qw(:,k - p - 1,:) = Qw(:,k - p,:);
        
        for l = 1 : p
            ind   = k - p + l;
            alpha = UQ(k + l) - X(j + 1);
            if abs(alpha) == 0
                Qw(:,ind - 1,:) = Qw(:,ind,:);
            else
                alpha = alpha / (UQ(k + l) - UP(i - p + l));
                Qw(:,ind - 1,:) = alpha * Qw(:,ind - 1,:) + (1 - alpha) * Qw(:,ind,:);
            end
        end
        
        UQ(k) = X(j + 1);
        k = k - 1;
    end
    %-Update
    Nurbs2D.cPoints{kk} = Qw;
    Nurbs2D.knots.U{kk} = transpose(UQ);
    Nurbs2D.number{kk}(1) = no2;
    [INC,IEN,nel,nnp,nen] = connectivity(Nurbs2D.order{kk},Nurbs2D.number{kk});
    Nurbs2D.INC{kk} = INC;
    Nurbs2D.IEN{kk} = IEN;
    Nurbs2D.nel{kk} = nel;
    Nurbs2D.nnp{kk} = nnp;
    Nurbs2D.nen{kk} = nen;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if vr == 0
    disp(' ');
elseif vr > 0
    X = linspace(0,1,vr+2);
    X = X(2:end-1);
    UP = Nurbs2D.knots.V{kk};
    Pw = Nurbs2D.cPoints{kk};
    mp = numel(Nurbs2D.knots.V{kk}) - 1;
    p  = Nurbs2D.order{kk}(2)-1;
    np = mp - p - 1;
    r  = numel(X);
    %
    
    no1 = numel(UP) + numel(X);
    UQ  = zeros(1, no1);
    
    no2 = size(Pw,2);
    no3 = size(Pw,3) + numel(X);
    Qw  = zeros(4, no2, no3);
    
    a = findspan(np,p,X(1),UP) - 1; %* -1
    b = findspan(np,p,X(end),UP);   %* -1
    
    for i = 0 : a - p
        Qw(:, :, i + 1) = Pw(:, :, i + 1);
    end
    
    for i = b - 1 : np
        Qw(:, :, i + r + 1) = Pw(:, :, i + 1);
    end
    
    for i = 0 : a
        UQ(i + 1) = UP(i + 1);
    end
    
    for i = b + p : mp
        UQ(i + r + 1) = UP(i + 1);
    end
    
    i = b + p - 1 + 1;
    k = b + p + r;
    
    for j = r - 1 : -1 : 0
        while X(j + 1) <= UP(i) && i - 1 > a
            Qw(:, :, k - p - 1) = Pw(:, :, i - p - 1);
            UQ(k) = UP(i);
            k = k - 1;
            i = i - 1;
        end
        
        Qw(:, :, k - p - 1) = Qw(:, :, k - p);
        
        for l = 1 : p
            ind   = k - p + l;
            alpha = UQ(k + l) - X(j + 1);
            if abs(alpha) == 0
                Qw(:, :, ind - 1) = Qw(:, :, ind);
            else
                alpha = alpha / (UQ(k + l) - UP(i - p + l));
                Qw(:, :, ind - 1) = alpha * Qw(:, :, ind - 1) + (1 - alpha) * Qw(:, :, ind);
            end
        end
        UQ(k) = X(j + 1);
        k = k - 1;
    end
    %-Update
    Nurbs2D.cPoints{kk} = Qw;
    Nurbs2D.knots.V{kk} = transpose(UQ);
    Nurbs2D.number{kk}(2) = no3;
    [INC,IEN,nel,nnp,nen] = connectivity(Nurbs2D.order{kk},Nurbs2D.number{kk});
    Nurbs2D.INC{kk} = INC;
    Nurbs2D.IEN{kk} = IEN;
    Nurbs2D.nel{kk} = nel;
    Nurbs2D.nnp{kk} = nnp;
    Nurbs2D.nen{kk} = nen;
end