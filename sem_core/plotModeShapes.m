function plotModeShapes(sem2D,modeNumPlot)
%
ldof = sem2D.local_dof;
%
for i = 1:modeNumPlot
    figure
    hold on
    
    if max(sem2D.uModes(1:3:end,i)) > -min(sem2D.uModes(1:3:end,i))
        modesign = 1;
        clim1 = max(sem2D.uModes(1:3:end,i));
    else
        modesign = -1;
        clim1 = -min(sem2D.uModes(1:3:end,i));
    end
    
    for el = 1:sem2D.nel
        %
        nconn = sem2D.conn(el,ldof:ldof:end)./ldof;
        %
        x_el = reshape(sem2D.nodes(nconn,1),sem2D.N,sem2D.N);
        y_el = reshape(sem2D.nodes(nconn,2),sem2D.N,sem2D.N);
        u_el = reshape(sem2D.uModes(sem2D.conn(el,1:3:end),i),sem2D.N,sem2D.N);
        %
        surf(x_el,y_el,modesign.*u_el)
    end
    %
    hold off
    axis equal
    axis off
    box on
    shading interp
    colormap jet
    clim([-clim1 clim1]);
    view(0,90)
    title(['Mode ' num2str(i)],'FontSize',12,'FontWeight','normal')
end