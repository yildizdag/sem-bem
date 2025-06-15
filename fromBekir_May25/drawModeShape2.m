
modenumber = 3;


indW = find((sem_mesh.ind_ALL(:,1) == 3));

figure
hold on

if max(eigVec(indW,modenumber)) > -min(eigVec(indW,modenumber))
    modesign = 1;
    clim1 = max(eigVec(indW,modenumber));
else
    modesign = -1;
    clim1 = -min(eigVec(indW,modenumber));
end

for di1 = 1:size(sem_mesh.elements,1)
    if sum(abs(sem_mesh.nodes(sem_mesh.elements(di1,:),3))) == 0
        
        pointsnow = sem_mesh.elementpoints(di1,:);
        pointsnow = pointsnow(pointsnow>0);
        
        xelm = zeros(sem_mesh.polynums(di1,2),sem_mesh.polynums(di1,1));        
        xelm(:) = round(sem_mesh.posn0(pointsnow,1),6);
        yelm = zeros(sem_mesh.polynums(di1,2),sem_mesh.polynums(di1,1));        
        yelm(:) = round(sem_mesh.posn0(pointsnow,2),6);
        
        modeshape = zeros(length(pointsnow),1);
        for di2 = 1:length(pointsnow)
            aa = find((sem_mesh.ind_ALL(:,3)==pointsnow(di2))&(sem_mesh.ind_ALL(:,1)==3));
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
colormap jet
caxis([-clim1 clim1])
view(0,90)
title(['Mode ' num2str(modenumber)],'FontSize',12,'FontWeight','normal')


