function  [xb, yb, dxbdxi, dybdxi, dxbdeta, dybdeta, fitx, fity] = ...
    Cross_section_Mapping (Mapping_Order, loc, xi, eta)
% -------------------------------------------------------------------------
% This code maps the points in a given cross-section into a square cross-section
% element (master element). The important point in mapping is that it
% should be checked whether a one-to-one mapping exists using the Jacobian
% matrix of the mapping. 
%
% ------------------------------------------------------------------------

fitx=zeros(5,5);
fity=zeros(5,5);

[psixi, psieta, constants, locA] = Mapping_Order_Phi(Mapping_Order);

% % figure(1)
% % % eval(['load loc',geometry]);
% % 
% % % Plotting the nodes only and numbering them
% % plot(loc(:,1),loc(:,2),'b.')
% % text(loc(:,1),loc(:,2),num2str(locA'),'FontSize',10)
% % axis equal

% Drawing the boundaries of each cross-section     
xisA=xi;
etasA=eta;
edge1x=zeros(length(xisA),1);
edge1y=zeros(length(xisA),1);
edge2x=zeros(length(xisA),1);
edge2y=zeros(length(xisA),1);
edge3x=zeros(length(etasA),1);
edge3y=zeros(length(etasA),1);
edge4x=zeros(length(etasA),1);
edge4y=zeros(length(etasA),1);

for i=1:length(xisA)
    for ii=1:25
            ii;
            edge1x(i,1)=edge1x(i,1)+constants(ii)*loc(ii,1)*polyval(psixi(ii,:),xisA(i))*polyval(psieta(ii,:),1);
            edge1y(i,1)=edge1y(i,1)+constants(ii)*loc(ii,2)*polyval(psixi(ii,:),xisA(i))*polyval(psieta(ii,:),1);

            edge2x(i,1)=edge2x(i,1)+constants(ii)*loc(ii,1)*polyval(psixi(ii,:),xisA(i))*polyval(psieta(ii,:),-1);
            edge2y(i,1)=edge2y(i,1)+constants(ii)*loc(ii,2)*polyval(psixi(ii,:),xisA(i))*polyval(psieta(ii,:),-1);
    end
end

for j=1:length(etasA)
    for ii=1:25
            ii;
            edge3x(j,1)=edge3x(j,1)+constants(ii)*loc(ii,1)*polyval(psixi(ii,:),1)*polyval(psieta(ii,:),etasA(j));
            edge3y(j,1)=edge3y(j,1)+constants(ii)*loc(ii,2)*polyval(psixi(ii,:),1)*polyval(psieta(ii,:),etasA(j));

            edge4x(j,1)=edge4x(j,1)+constants(ii)*loc(ii,1)*polyval(psixi(ii,:),-1)*polyval(psieta(ii,:),etasA(j));
            edge4y(j,1)=edge4y(j,1)+constants(ii)*loc(ii,2)*polyval(psixi(ii,:),-1)*polyval(psieta(ii,:),etasA(j));
    end
end
% % % hold on
% % % plot(edge1x,edge1y,'b',edge2x,edge2y,'r',edge3x,edge3y,'k',edge4x,edge4y,'y')
% % % hold off
% Mapping cross-section information (polynomial coefficients) into a 3D tensor 
for ii=1:25
    for jj=1:5
        for kk=1:5
            fitx(jj,kk)=fitx(jj,kk)+loc(ii,1)*constants(ii)*psixi(ii,jj)*psieta(ii,kk);
            fity(jj,kk)=fity(jj,kk)+loc(ii,2)*constants(ii)*psixi(ii,jj)*psieta(ii,kk);
        end
    end 
end
    


% Finding the nodes at each cross-section and constructing them in a tensor form
xb = zeros(length(xi),length(eta));
yb = zeros(length(xi),length(eta));
zb = zeros(length(xi),length(eta));

% % figure(2)
mapPointsX=zeros(length(xisA),length(etasA));
mapPointsY=zeros(length(xisA),length(etasA));
mapPointsZ=zeros(length(xisA),length(etasA));
for i=1:length(xisA)
    for j=1:length(etasA)
        for ii=1:5
            for jj=1:5

                mapPointsX(i,j) = mapPointsX(i,j)+fitx(ii,jj)*xisA(i)^(5-ii)*etasA(j)^(5-jj);
                mapPointsY(i,j) = mapPointsY(i,j)+fity(ii,jj)*xisA(i)^(5-ii)*etasA(j)^(5-jj);
                mapPointsZ(i,j) = 0;

                xb(i,j)=xb(i,j)+fitx(ii,jj)*xisA(i)^(5-ii)*etasA(j)^(5-jj);
                yb(i,j)=yb(i,j)+fity(ii,jj)*xisA(i)^(5-ii)*etasA(j)^(5-jj);
            end
        end
        zb(i,j)=0;
    end
end
% % plot3(mapPointsX, mapPointsY, mapPointsZ,'r.')
% % hold on
% % axis equal



% Calculating the derivatives to transform into the mapped domain
dxbdxi=zeros(length(xisA),length(etasA));
for i=1:length(xisA)
    for j=1:length(etasA)
        for ii=1:5
            for jj=1:5
                if ii==5
                    continue
                end
                dxbdxi(i,j)=dxbdxi(i,j)+fitx(ii,jj)*(5-ii)*xisA(i)^(4-ii)*etasA(j)^(5-jj);
            end
        end
    end
end

dybdxi=zeros(length(xisA),length(etasA));
for i=1:length(xisA)
    for j=1:length(etasA)
        for ii=1:5
            for jj=1:5
                if ii==5
                    continue
                end
                dybdxi(i,j)=dybdxi(i,j)+fity(ii,jj)*(5-ii)*xisA(i)^(4-ii)*etasA(j)^(5-jj);
            end
        end
    end
end

dxbdeta=zeros(length(xisA),length(etasA));
for i=1:length(xisA)
    for j=1:length(etasA)
        for ii=1:5
            for jj=1:5
                if jj==5
                    continue
                end
                dxbdeta(i,j)=dxbdeta(i,j)+fitx(ii,jj)*xisA(i)^(5-ii)*(5-jj)*etasA(j)^(4-jj);
            end
        end
    end
end


dybdeta=zeros(length(xisA),length(etasA));
for i=1:length(xisA)
    for j=1:length(etasA)
        for ii=1:5
            for jj=1:5
                if jj==5
                    continue
                end
                dybdeta(i,j)=dybdeta(i,j)+fity(ii,jj)*xisA(i)^(5-ii)*(5-jj)*etasA(j)^(4-jj);
            end
        end
    end
end


