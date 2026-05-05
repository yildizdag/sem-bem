function [Nurbs2D_plate,Nurbs2D_stiff] = update_baseline(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp)
%
Nurbs2D_plate.nodesDeformed = Nurbs2D_plate.nodes; %copy original control points
Nurbs2D_stiff.nodesDeformed = Nurbs2D_stiff.nodes; %copy original control points
%
for j = 1:size(pconn,1)
    p1 = pconn(j,1); e1 = pconn(j,2);
    p2 = pconn(j,3); e2 = pconn(j,4);
    s  = pconn(j,5);
    numMCP = size(Nurbs2D_plate.movingCP{p1,e1},2);
    indMCP1 = sort(Nurbs2D_plate.movingCP{p1,e1});
    indMCP2 = sort(Nurbs2D_plate.movingCP{p2,e2});
    dispx = dcp; %randomly generated displacement along x: -xrange < dispx < xrange
    dispy = 0*(2.*rand(1,numMCP)-1); %randomly generated displacement along y: -yrange < dispy < yrange
    for k = 1:numMCP
        indStiff = find(ismember(Nurbs2D_stiff.nodes{s}(:,1),Nurbs2D_plate.nodes{p1}(indMCP1(k),1))...
                      & ismember(Nurbs2D_stiff.nodes{s}(:,2),Nurbs2D_plate.nodes{p1}(indMCP1(k),2)));
        Nurbs2D_stiff.nodesDeformed{s}(indStiff,1) = Nurbs2D_stiff.nodes{s}(indStiff,1) + transpose(dispx(k));
        Nurbs2D_stiff.nodesDeformed{s}(indStiff,2) = Nurbs2D_stiff.nodes{s}(indStiff,2) + transpose(dispy(k));
    end
    Nurbs2D_stiff.cPoints{s} = reshape(transpose(Nurbs2D_stiff.nodesDeformed{s}),4,Nurbs2D_stiff.number{s}(1,1),Nurbs2D_stiff.number{s}(1,2));
    Nurbs2D_plate.nodesDeformed{p1}(indMCP1,1) = Nurbs2D_plate.nodes{p1}(indMCP1,1) + transpose(dispx);
    Nurbs2D_plate.nodesDeformed{p1}(indMCP1,2) = Nurbs2D_plate.nodes{p1}(indMCP1,2) + transpose(dispy);
    Nurbs2D_plate.cPoints{p1} = reshape(transpose(Nurbs2D_plate.nodesDeformed{p1}),4,Nurbs2D_plate.number{p1}(1,1),Nurbs2D_plate.number{p1}(1,2));
    Nurbs2D_plate.nodesDeformed{p2}(indMCP2,1) = Nurbs2D_plate.nodes{p2}(indMCP2,1) + transpose(dispx);
    Nurbs2D_plate.nodesDeformed{p2}(indMCP2,2) = Nurbs2D_plate.nodes{p2}(indMCP2,2) + transpose(dispy);
    Nurbs2D_plate.cPoints{p2} = reshape(transpose(Nurbs2D_plate.nodesDeformed{p2}),4,Nurbs2D_plate.number{p2}(1,1),Nurbs2D_plate.number{p2}(1,2));
end

