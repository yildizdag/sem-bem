
modenumber = 2;


indW = find((indA(:,1) == 3));

figure
hold on

if max(eigVec(indW,modenumber)) > -min(eigVec(indW,modenumber))
    modesign = 1;
    clim1 = max(eigVec(indW,modenumber));
else
    modesign = -1;
    clim1 = -min(eigVec(indW,modenumber));
end

for di1 = 1:size(elements,1)
    if sum(abs(nodes(elements(di1,:),3))) == 0
        
        pointsnow = elementpoints(di1,:);
        pointsnow = pointsnow(pointsnow>0);
        
        xelm = zeros(polynum_eta(di1),polynum_xi(di1));        
        xelm(:) = round(posn0(pointsnow,1),6);
        yelm = zeros(polynum_eta(di1),polynum_xi(di1));        
        yelm(:) = round(posn0(pointsnow,2),6);
        
        modeshape = zeros(length(pointsnow),1);
        for di2 = 1:length(pointsnow)
            aa = find((indB==pointsnow(di2))&(indA==3));
            if ~isempty(aa)
                modeshape(di2) = eigVec(aa,modenumber);
            end            
        end
        
        modeshape = reshape(modeshape,size(xelm));
        
        surf(xelm,yelm,modeshape*modesign)
    end
end

axis equal
axis off
box on
shading interp
colormap summer
caxis([-clim1 clim1])
view(0,90)
title(['Mode ' num2str(modenumber)],'FontSize',12,'FontWeight','normal')


