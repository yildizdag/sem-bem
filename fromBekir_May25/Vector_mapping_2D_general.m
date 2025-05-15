function [VD, QDx, QDy] = Vector_mapping_2D_general (polynum_xi,...
    polynum_eta, xi, eta, V_xi, V_eta, Q1_xi, Q1_eta, JAC, ...
    dxidx, detadx, dxidy, detady, space_xi, space_eta)

% ------------------------------------------------------------------------

spn_xi = polynum_xi;
spn_eta = polynum_eta;

spn_xy = spn_xi*spn_eta;

% ---- Expanding/Mapping of Inner Product Matrix ----
% VD = Fxy_mapping (spn_xi, spn_eta, V_xi, V_eta);
% VD = InnerProduct_3functions(space_xi,space_eta,JAC);

VD = zeros(spn_xy,spn_xy);
for j=1:spn_xi
    for k=1:spn_eta        
        for j2=1:spn_xi
            for k2=1:spn_eta
                count = (j-1)*spn_eta+k;
                count2 = (j2-1)*spn_eta+k2;
                VD(count2,count) = VD(count2,count) + JAC(j,k)*V_xi(j2,j)*V_eta(k2,k);
            end
        end
    end
end


% ---- Mapping (Expanding) of Derivative Matrices ----
QDxi_dxidx=zeros(spn_xy,spn_xy);
for i = 1:length(xi)
    for j = 1:length(eta)
        for i2 = 1:length(xi)
            count = (i-1)*spn_eta+j;
            count2 = (i2-1)*spn_eta+j;
            QDxi_dxidx(count2,count) = QDxi_dxidx(count2,count)+...
                dxidx(i2,j)*Q1_xi(i2,i);
        end
    end
end

QDeta_detadx = zeros(spn_xy,spn_xy);
for i = 1:length(xi)
    for j = 1:length(eta)
        for j2 = 1:length(eta)
            count = (i-1)*spn_eta+j;
            count2 = (i-1)*spn_eta+j2;
            QDeta_detadx(count2,count) = QDeta_detadx(count2,count)+...
                detadx(i,j2)*Q1_eta(j2,j);
        end
    end
end

QDxi_dxidy = zeros(spn_xy,spn_xy);
for i = 1:length(xi)
    for j = 1:length(eta)
        for i2 = 1:length(xi)
            count = (i-1)*spn_eta+j;
            count2 = (i2-1)*spn_eta+j;
            QDxi_dxidy(count2,count) = QDxi_dxidy(count2,count)+...
                Q1_xi(i2,i)*dxidy(i2,j);
        end
    end
end

QDeta_detady = zeros(spn_xy,spn_xy);
for i = 1:length(xi)
    for j = 1:length(eta)
        for j2 = 1:length(eta)
            count = (i-1)*spn_eta+j;
            count2 = (i-1)*spn_eta+j2;
            QDeta_detady(count2,count) = QDeta_detady(count2,count)+...
                Q1_eta(j2,j)*detady(i,j2);
        end
    end
end

QDx = QDxi_dxidx + QDeta_detadx;
QDy = QDxi_dxidy + QDeta_detady;


end








