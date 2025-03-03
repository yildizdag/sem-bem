function iga2DmeshPlotNURBS(Nurbs2D)
%---------------------------------------------------------------
% Plots 2-D NURBS Patches
%---------------------------------------------------------------
for kk = 1:Nurbs2D.numDryPatch
    %Sampling:
    N1 = Nurbs2D.knots.U{kk}(end);
    N2 = Nurbs2D.knots.V{kk}(end);
    sampleU = ceil(N1/0.05); uu = linspace(Nurbs2D.knots.U{kk}(Nurbs2D.order{kk}(1)),N1,sampleU);
    sampleV = ceil(N2/0.05); vv = linspace(Nurbs2D.knots.V{kk}(Nurbs2D.order{kk}(2)),N2,sampleV);
    %Storing:
    X0 = zeros(length(vv),length(uu));
    Y0 = zeros(length(vv),length(uu));
    Z0 = zeros(length(vv),length(uu));
    for j = 1:sampleV
        for r = 1:sampleU
            iu = findspan(Nurbs2D.number{kk}(1),Nurbs2D.order{kk}(1)-1,uu(r),Nurbs2D.knots.U{kk});
            iv = findspan(Nurbs2D.number{kk}(2),Nurbs2D.order{kk}(2)-1,vv(j),Nurbs2D.knots.V{kk});
            dNu = dersbasisfuns(iu,uu(r),Nurbs2D.order{kk}(1)-1,0,Nurbs2D.knots.U{kk});
            dNv = dersbasisfuns(iv,vv(j),Nurbs2D.order{kk}(2)-1,0,Nurbs2D.knots.V{kk});
            CP = Nurbs2D.cPoints{kk}(:,iu-Nurbs2D.order{kk}(1)+1:iu, iv-Nurbs2D.order{kk}(2)+1:iv);
            Sw = zeros(4,1);
            for i = 1:4
                Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{kk}(1),Nurbs2D.order{kk}(2))*dNv(1,:)';
            end
            S = Sw(1:3,:) / Sw(4);
            X0(j,r) = S(1); Y0(j,r) = S(2); Z0(j,r) = S(3);
        end
    end
    hold on
    surf(X0,Y0,Z0,'EdgeLighting','none','FaceColor','[.8,1,.8]','EdgeColor','none','LineStyle','-');
    hold off
end

hold on
for kk = 1:Nurbs2D.numDryPatch
    xParametric = Nurbs2D.knots.U{kk}(Nurbs2D.order{kk}(1):end-Nurbs2D.order{kk}(1)+1);
    yParametric = Nurbs2D.knots.V{kk}(Nurbs2D.order{kk}(2):end-Nurbs2D.order{kk}(2)+1);
    N1 = Nurbs2D.knots.U{kk}(end);
    N2 = Nurbs2D.knots.V{kk}(end);
    sampleU = ceil(N1/0.05); uu = linspace(Nurbs2D.knots.U{kk}(Nurbs2D.order{kk}(1)),N1,sampleU);
    sampleV = ceil(N2/0.05); vv = linspace(Nurbs2D.knots.V{kk}(Nurbs2D.order{kk}(2)),N2,sampleV);
    for ii=1:length(yParametric)
        X=zeros(1,sampleU); Y=zeros(1,sampleU); Z=zeros(1,sampleU);
        for jj=1:sampleU
            iu = findspan(Nurbs2D.number{kk}(1),Nurbs2D.order{kk}(1)-1,uu(jj),Nurbs2D.knots.U{kk});
            iv = findspan(Nurbs2D.number{kk}(2),Nurbs2D.order{kk}(2)-1,yParametric(ii),Nurbs2D.knots.V{kk});
            dNu = dersbasisfuns(iu,uu(jj),Nurbs2D.order{kk}(1)-1,0,Nurbs2D.knots.U{kk});
            dNv = dersbasisfuns(iv,yParametric(ii),Nurbs2D.order{kk}(2)-1,0,Nurbs2D.knots.V{kk});
            CP = Nurbs2D.cPoints{kk}(:,iu-Nurbs2D.order{kk}(1)+1:iu, iv-Nurbs2D.order{kk}(2)+1:iv);
            Sw = zeros(4,1);
            for i = 1:4
                Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{kk}(1),Nurbs2D.order{kk}(2))*dNv(1,:)';
            end
            S = Sw(1:3,:) / Sw(4);
            X(jj) = S(1); Y(jj) = S(2); Z(jj) = S(3);
        end
        plot3(X,Y,Z,'k','LineWidth',3);
    end
    for ii=1:length(xParametric)
        X=zeros(1,sampleV); Y=zeros(1,sampleV); Z=zeros(1,sampleV);
        for jj=1:sampleV
            iu = findspan(Nurbs2D.number{kk}(1),Nurbs2D.order{kk}(1)-1,xParametric(ii),Nurbs2D.knots.U{kk});
            iv = findspan(Nurbs2D.number{kk}(2),Nurbs2D.order{kk}(2)-1,vv(jj),Nurbs2D.knots.V{kk});
            dNu = dersbasisfuns(iu,xParametric(ii),Nurbs2D.order{kk}(1)-1,0,Nurbs2D.knots.U{kk});
            dNv = dersbasisfuns(iv,vv(jj),Nurbs2D.order{kk}(2)-1,0,Nurbs2D.knots.V{kk});
            CP = Nurbs2D.cPoints{kk}(:,iu-Nurbs2D.order{kk}(1)+1:iu, iv-Nurbs2D.order{kk}(2)+1:iv);
            Sw = zeros(4,1);
            for i = 1:4
                Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{kk}(1),Nurbs2D.order{kk}(2))*dNv(1,:)';
            end
            S = Sw(1:3,:) / Sw(4);
            X(jj) = S(1); Y(jj) = S(2); Z(jj) = S(3);
        end
        plot3(X,Y,Z,'k','LineWidth',3);
    end
end
hold off
axis equal