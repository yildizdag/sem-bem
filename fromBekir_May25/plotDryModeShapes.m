function plotDryModeShapes(sem_mesh,numDryModes)
%
% Dry Mode Shapes (SEM Results)
% Generates modal shapes using ONLY sampling points data
% Further improvement is required to interpolate within each element
%
for i = 1:numDryModes
    %
    indW = find((sem_mesh.ind_ALL(:,1) == 3));
    %
    figure
    hold on
    %
    if max(sem_mesh.eigVec(indW,i)) > -min(sem_mesh.eigVec(indW,i))
        modesign = 1;
        clim1 = max(sem_mesh.eigVec(indW,i));
    else
        modesign = -1;
        clim1 = -min(sem_mesh.eigVec(indW,i));
    end
    %
    for di1 = 1:size(sem_mesh.elements,1)
        if sum(abs(sem_mesh.nodes(sem_mesh.elements(di1,:),3))) == 0
            %
            pointsnow = sem_mesh.elementpoints(di1,:);
            pointsnow = pointsnow(pointsnow>0);
            %
            xelm = zeros(sem_mesh.polynums(di1,2),sem_mesh.polynums(di1,1));
            xelm(:) = round(sem_mesh.posn0(pointsnow,1),6);
            yelm = zeros(sem_mesh.polynums(di1,2),sem_mesh.polynums(di1,1));
            yelm(:) = round(sem_mesh.posn0(pointsnow,2),6);
            %
            modeshape = zeros(length(pointsnow),1);
            for di2 = 1:length(pointsnow)
                aa = find((sem_mesh.ind_ALL(:,3)==pointsnow(di2))&(sem_mesh.ind_ALL(:,1)==3));
                if ~isempty(aa)
                    modeshape(di2) = sem_mesh.eigVec(aa,i);
                end
            end
            %
            modeshape = reshape(modeshape,size(xelm));
            %
            surf(xelm,yelm,modeshape*modesign)
        end
    end
    %
    axis equal
    axis off
    box on
    shading interp
    colormap jet
    caxis([-clim1 clim1])
    view(0,90)
    title(['Mode ' num2str(i)],'FontSize',12,'FontWeight','normal')
end


