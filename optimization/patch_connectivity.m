function [Nurbs2D_plate,Nurbs2D_stiff,pconn] = patch_connectivity(Nurbs2D_plate,Nurbs2D_stiff)
pconn = [];
ccount = 1;
outer = 0;
for p1 = 1:Nurbs2D_plate.numpatch
    p1_edge = [Nurbs2D_plate.edges.coefsNo{p1,1}(1),Nurbs2D_plate.edges.coefsNo{p1,1}(end);
               Nurbs2D_plate.edges.coefsNo{p1,2}(1),Nurbs2D_plate.edges.coefsNo{p1,2}(end);
               Nurbs2D_plate.edges.coefsNo{p1,3}(1),Nurbs2D_plate.edges.coefsNo{p1,3}(end);
               Nurbs2D_plate.edges.coefsNo{p1,4}(1),Nurbs2D_plate.edges.coefsNo{p1,4}(end)];
    p1_cp = zeros(2,3,4);
    p1_cp(1,:,1) = Nurbs2D_plate.nodes{p1}(p1_edge(1,1),1:3); p1_cp(2,:,1) = Nurbs2D_plate.nodes{p1}(p1_edge(1,2),1:3);
    p1_cp(1,:,2) = Nurbs2D_plate.nodes{p1}(p1_edge(2,1),1:3); p1_cp(2,:,2) = Nurbs2D_plate.nodes{p1}(p1_edge(2,2),1:3);
    p1_cp(1,:,3) = Nurbs2D_plate.nodes{p1}(p1_edge(3,1),1:3); p1_cp(2,:,3) = Nurbs2D_plate.nodes{p1}(p1_edge(3,2),1:3);
    p1_cp(1,:,4) = Nurbs2D_plate.nodes{p1}(p1_edge(4,1),1:3); p1_cp(2,:,4) = Nurbs2D_plate.nodes{p1}(p1_edge(4,2),1:3);
    for p2 = (p1+1):Nurbs2D_plate.numpatch
        p2_edge = [Nurbs2D_plate.edges.coefsNo{p2,1}(1),Nurbs2D_plate.edges.coefsNo{p2,1}(end);
                   Nurbs2D_plate.edges.coefsNo{p2,2}(1),Nurbs2D_plate.edges.coefsNo{p2,2}(end);
                   Nurbs2D_plate.edges.coefsNo{p2,3}(1),Nurbs2D_plate.edges.coefsNo{p2,3}(end);
                   Nurbs2D_plate.edges.coefsNo{p2,4}(1),Nurbs2D_plate.edges.coefsNo{p2,4}(end)];
        p2_cp = zeros(2,3,4);
        p2_cp(1,:,1) = Nurbs2D_plate.nodes{p2}(p2_edge(1,1),1:3); p2_cp(2,:,1) = Nurbs2D_plate.nodes{p2}(p2_edge(1,2),1:3);
        p2_cp(1,:,2) = Nurbs2D_plate.nodes{p2}(p2_edge(2,1),1:3); p2_cp(2,:,2) = Nurbs2D_plate.nodes{p2}(p2_edge(2,2),1:3);
        p2_cp(1,:,3) = Nurbs2D_plate.nodes{p2}(p2_edge(3,1),1:3); p2_cp(2,:,3) = Nurbs2D_plate.nodes{p2}(p2_edge(3,2),1:3);
        p2_cp(1,:,4) = Nurbs2D_plate.nodes{p2}(p2_edge(4,1),1:3); p2_cp(2,:,4) = Nurbs2D_plate.nodes{p2}(p2_edge(4,2),1:3);
        for i = 1:4
            CP1 = p1_cp(:,:,i);
            for j = 1:4
                CP2 = p2_cp(:,:,j);
                checkCP = ismembertol(CP1,CP2,'ByRows',true);
                if checkCP(1) && checkCP(2)
                    pconn(ccount,:) = [p1 i p2 j 0 0];
                    if outer == 1
                        Nurbs2D_plate.movingCP{p1,i} = Nurbs2D_plate.edges.coefsNo{p1,i}(1:end);
                        Nurbs2D_plate.movingCP{p2,j} = Nurbs2D_plate.edges.coefsNo{p2,j}(1:end);
                    else
                        Nurbs2D_plate.movingCP{p1,i} = Nurbs2D_plate.edges.coefsNo{p1,i}(2:end-1);
                        Nurbs2D_plate.movingCP{p2,j} = Nurbs2D_plate.edges.coefsNo{p2,j}(2:end-1);
                    end
                    ccount = ccount + 1;
                end
            end
        end
    end
end
% ----------------------------------
% Patch Connectivity for Stiffeners
% ----------------------------------
for s = 1:Nurbs2D_stiff.numpatch
    s_edge = [Nurbs2D_stiff.edges.coefsNo{s,1}(1),Nurbs2D_stiff.edges.coefsNo{s,1}(end);
              Nurbs2D_stiff.edges.coefsNo{s,2}(1),Nurbs2D_stiff.edges.coefsNo{s,2}(end);
              Nurbs2D_stiff.edges.coefsNo{s,3}(1),Nurbs2D_stiff.edges.coefsNo{s,3}(end);
              Nurbs2D_stiff.edges.coefsNo{s,4}(1),Nurbs2D_stiff.edges.coefsNo{s,4}(end)];
    s_cp = zeros(2,3,4);
    s_cp(1,:,1) = Nurbs2D_stiff.nodes{s}(s_edge(1,1),1:3); s_cp(2,:,1) = Nurbs2D_stiff.nodes{s}(s_edge(1,2),1:3);
    s_cp(1,:,2) = Nurbs2D_stiff.nodes{s}(s_edge(2,1),1:3); s_cp(2,:,2) = Nurbs2D_stiff.nodes{s}(s_edge(2,2),1:3);
    s_cp(1,:,3) = Nurbs2D_stiff.nodes{s}(s_edge(3,1),1:3); s_cp(2,:,3) = Nurbs2D_stiff.nodes{s}(s_edge(3,2),1:3);
    s_cp(1,:,4) = Nurbs2D_stiff.nodes{s}(s_edge(4,1),1:3); s_cp(2,:,4) = Nurbs2D_stiff.nodes{s}(s_edge(4,2),1:3);
    for p = 1:size(pconn,1)
        p_edge = [Nurbs2D_plate.edges.coefsNo{pconn(p,1),pconn(p,2)}(1),Nurbs2D_plate.edges.coefsNo{pconn(p,1),pconn(p,2)}(end)];
        p_cp = [Nurbs2D_plate.nodes{pconn(p,1)}(p_edge(1),1:3); Nurbs2D_plate.nodes{pconn(p,1)}(p_edge(2),1:3)];
        for i = 1:4
            checkCP = ismembertol(s_cp(:,:,i),p_cp,'ByRows',true);
            if checkCP(1) && checkCP(2)
                pconn(p,5:6) = [s i];
                if outer == 1
                    Nurbs2D_stiff.movingCP{s,i} = Nurbs2D_stiff.edges.coefsNo{s,i}(1:end);
                else
                    Nurbs2D_stiff.movingCP{s,i} = Nurbs2D_stiff.edges.coefsNo{s,i}(2:end-1);
                end
            end
        end
    end
end