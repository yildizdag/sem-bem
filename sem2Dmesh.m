function sem2Dmesh(FileName,semPatch,bemPatch,np_u,np_v,plotNURBS,plotSEM)
% NURBS-Enhanced coarse quad mesh generator
% Process provided NURBS data
numSEMpatch = size(semPatch,2);
numBEMpatch = size(bemPatch,2);
%------------------------------------------------------------------
% Create 2-D Nurbs Structure (reads FileName)
Nurbs2D = iga2DmeshDry(FileName,numSEMpatch,1);
%------------------------------------------------------------------
% Plot Imported 2-D NURBS Structure
if plotNURBS == 1
    figure;
    iga2DmeshPlotNURBS(Nurbs2D);
    axis off
end
%------------------------------------------------------------------
% Points for Spectral Element Method (e.g. 5 x 5, 3 x 3, etc.)
tot_el = 0; %Total num of elements
for k = 1:Nurbs2D.numDryPatch
    tot_el = tot_el + Nurbs2D.nel{k};
end
elData = zeros(np_u*np_v,3,tot_el);
nodeData = zeros(np_u*np_v*tot_el,3);
count_el = 1;
count_node = 1;
if plotSEM == 1
    figure;
    iga2DmeshPlotNURBS(Nurbs2D);
    for k = 1:Nurbs2D.numDryPatch
        for el = 1:Nurbs2D.nel{k}
            iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
            iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
            u_sample = linspace(Nurbs2D.knots.U{k}(iu),Nurbs2D.knots.U{k}(iu+1),np_u);
            v_sample = linspace(Nurbs2D.knots.V{k}(iv),Nurbs2D.knots.V{k}(iv+1),np_v);
            x_sample = zeros(np_u*np_v,1);
            y_sample = zeros(np_u*np_v,1);
            z_sample = zeros(np_u*np_v,1);
            count = 1;
            for j = 1:np_u
                for r = 1:np_v
                    dNu = dersbasisfuns(iu,u_sample(j),Nurbs2D.order{k}(1)-1,0,Nurbs2D.knots.U{k});
                    dNv = dersbasisfuns(iv,v_sample(r),Nurbs2D.order{k}(2)-1,0,Nurbs2D.knots.V{k});
                    CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                    Sw = zeros(4,1);
                    for i = 1:4
                        Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{k}(1),Nurbs2D.order{k}(2))*dNv(1,:)';
                    end
                    S = Sw(1:3,:) / Sw(4);
                    x_sample(count) = S(1);
                    y_sample(count) = S(2);
                    z_sample(count) = S(3);
                    elData(count,:,count_el) = [S(1) S(2) S(3)];
                    nodeData(count_node,:) = [S(1), S(2), S(3)];
                    count = count+1;
                    count_node = count_node+1;
                end
            end
            hold on
            scatter3(x_sample,y_sample,z_sample,80,'b','filled');
            hold off
            count_el = count_el+1;
        end
    end
    axis off
end
%---------------------------------------------------------------
% Patch elementsectivity
% Nodal Coordinates (nodes)
% elementsectivity (elements)
%---------------------------------------------------------------
TOL = 0.001; %---> Check!
nodes = uniquetol(nodeData,TOL,'ByRows',true);
elements = zeros(tot_el,np_u*np_v);
for i = 1:tot_el
    for j = 1:np_u*np_v
        node_id = find(ismembertol(nodes, elData(j,:,i),TOL,'ByRows',true));
        elements(i,j) = node_id;
    end
end
save('nodes.mat','nodes');
save('elements.mat','elements');
%---------------------------------------------------------------
%---------------------------------------------------------------
% BEM MESH
Nurbs2D = iga2DmeshWet(FileName,numSEMpatch,numBEMpatch,1,Nurbs2D);
%---------------------------------------------------------------
tot_el = 0; %Total num of elements
for k = Nurbs2D.numDryPatch+1:Nurbs2D.numpatch
    tot_el = tot_el + Nurbs2D.nel{k};
end
pbem = 2;
elData = zeros(pbem*pbem,3,tot_el);
nodeData = zeros(pbem*pbem*tot_el,3);
count_el = 1;
count_node = 1;
if plotSEM == 1
    figure;
    iga2DmeshPlotNURBS(Nurbs2D);
    for k = Nurbs2D.numDryPatch+1:Nurbs2D.numpatch
        for el = 1:Nurbs2D.nel{k}
            iu = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),1);   
            iv = Nurbs2D.INC{k}(Nurbs2D.IEN{k}(1,el),2);
            u_sample = linspace(Nurbs2D.knots.U{k}(iu),Nurbs2D.knots.U{k}(iu+1),pbem);
            v_sample = linspace(Nurbs2D.knots.V{k}(iv),Nurbs2D.knots.V{k}(iv+1),pbem);
            x_sample = zeros(pbem*pbem,1);
            y_sample = zeros(pbem*pbem,1);
            z_sample = zeros(pbem*pbem,1);
            count = 1;
            for j = 1:pbem
                for r = 1:pbem
                    dNu = dersbasisfuns(iu,u_sample(j),Nurbs2D.order{k}(1)-1,0,Nurbs2D.knots.U{k});
                    dNv = dersbasisfuns(iv,v_sample(r),Nurbs2D.order{k}(2)-1,0,Nurbs2D.knots.V{k});
                    CP = Nurbs2D.cPoints{k}(:,iu-Nurbs2D.order{k}(1)+1:iu, iv-Nurbs2D.order{k}(2)+1:iv);
                    Sw = zeros(4,1);
                    for i = 1:4
                        Sw(i,1) = dNu(1,:)*reshape(CP(i,:,:),Nurbs2D.order{k}(1),Nurbs2D.order{k}(2))*dNv(1,:)';
                    end
                    S = Sw(1:3,:) / Sw(4);
                    x_sample(count) = S(1);
                    y_sample(count) = S(2);
                    z_sample(count) = S(3);
                    elData(count,:,count_el) = [S(1) S(2) S(3)];
                    nodeData(count_node,:) = [S(1), S(2), S(3)];
                    count = count+1;
                    count_node = count_node+1;
                end
            end
            hold on
            scatter3(x_sample,y_sample,z_sample,80,'b','filled');
            hold off
            count_el = count_el+1;
        end
    end
    axis off
end
%---------------------------------------------------------------
% Patch elementsectivity
% Nodal Coordinates (nodes)
% elementsectivity (elements)
%---------------------------------------------------------------
TOL = 0.001; %---> Check!
nodesBEM = uniquetol(nodeData,TOL,'ByRows',true);
elementsBEM = zeros(tot_el,pbem*pbem);
for i = 1:tot_el
    for j = 1:pbem*pbem
        node_id = find(ismembertol(nodesBEM, elData(j,:,i),TOL,'ByRows',true));
        elementsBEM(i,j) = node_id;
    end
end
save('nodesBEM.mat','nodesBEM');
save('elementsBEM.mat','elementsBEM');