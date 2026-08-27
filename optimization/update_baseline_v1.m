function [Nurbs2D_plate,Nurbs2D_stiff] = update_baseline_v1(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp)
%
% Minimal modification of the original update_baseline.m.
% dcp contains one x-displacement for each moving control-point group.
% The first and second rows of movingCP are updated separately to avoid
% incompatible array-size operations.
%

Nurbs2D_plate.nodesDeformed = Nurbs2D_plate.nodes;
Nurbs2D_stiff.nodesDeformed = Nurbs2D_stiff.nodes;

for j = 1:size(pconn,1)

    p1 = pconn(j,1); e1 = pconn(j,2);
    p2 = pconn(j,3); e2 = pconn(j,4);
    s  = pconn(j,5);

    indMCP1 = Nurbs2D_plate.movingCP{p1,e1};
    indMCP2 = Nurbs2D_plate.movingCP{p2,e2};

    numMCP = size(indMCP1,2);

    dispx = dcp(:).';
    dispy = zeros(1,numMCP);

    if numel(dispx) ~= numMCP
        error('update_baseline_v1:InvalidDesignVector', ...
              'dcp must contain exactly %d values.',numMCP);
    end

    for k = 1:numMCP

        % Interface control points
        cp1_interface = indMCP1(1,k);
        cp2_interface = indMCP2(1,k);

        % Adjacent control points, when present
        cp1_adjacent = [];
        cp2_adjacent = [];

        if size(indMCP1,1) >= 2
            cp1_adjacent = indMCP1(2,k);
        end

        if size(indMCP2,1) >= 2
            cp2_adjacent = indMCP2(2,k);
        end

        % Find stiffener control points from the interface point only
        indStiff = find( ...
            ismember(Nurbs2D_stiff.nodes{s}(:,1), ...
                     Nurbs2D_plate.nodes{p1}(cp1_interface,1)) & ...
            ismember(Nurbs2D_stiff.nodes{s}(:,2), ...
                     Nurbs2D_plate.nodes{p1}(cp1_interface,2)));

        % Update stiffener points
        for q = 1:numel(indStiff)
            id = indStiff(q);

            Nurbs2D_stiff.nodesDeformed{s}(id,1) = ...
                Nurbs2D_stiff.nodes{s}(id,1) + dispx(k);

            Nurbs2D_stiff.nodesDeformed{s}(id,2) = ...
                Nurbs2D_stiff.nodes{s}(id,2) + dispy(k);
        end

        % Update interface CP of plate patch 1
        Nurbs2D_plate.nodesDeformed{p1}(cp1_interface,1) = ...
            Nurbs2D_plate.nodes{p1}(cp1_interface,1) + dispx(k);

        Nurbs2D_plate.nodesDeformed{p1}(cp1_interface,2) = ...
            Nurbs2D_plate.nodes{p1}(cp1_interface,2) + dispy(k);

        % Update adjacent CP of plate patch 1
        if ~isempty(cp1_adjacent) && cp1_adjacent > 0
            Nurbs2D_plate.nodesDeformed{p1}(cp1_adjacent,1) = ...
                Nurbs2D_plate.nodes{p1}(cp1_adjacent,1) + dispx(k);

            Nurbs2D_plate.nodesDeformed{p1}(cp1_adjacent,2) = ...
                Nurbs2D_plate.nodes{p1}(cp1_adjacent,2) + dispy(k);
        end

        % Update interface CP of plate patch 2
        Nurbs2D_plate.nodesDeformed{p2}(cp2_interface,1) = ...
            Nurbs2D_plate.nodes{p2}(cp2_interface,1) + dispx(k);

        Nurbs2D_plate.nodesDeformed{p2}(cp2_interface,2) = ...
            Nurbs2D_plate.nodes{p2}(cp2_interface,2) + dispy(k);

        % Update adjacent CP of plate patch 2
        if ~isempty(cp2_adjacent) && cp2_adjacent > 0
            Nurbs2D_plate.nodesDeformed{p2}(cp2_adjacent,1) = ...
                Nurbs2D_plate.nodes{p2}(cp2_adjacent,1) + dispx(k);

            Nurbs2D_plate.nodesDeformed{p2}(cp2_adjacent,2) = ...
                Nurbs2D_plate.nodes{p2}(cp2_adjacent,2) + dispy(k);
        end
    end

    Nurbs2D_stiff.cPoints{s} = reshape( ...
        Nurbs2D_stiff.nodesDeformed{s}.', ...
        4,Nurbs2D_stiff.number{s}(1),Nurbs2D_stiff.number{s}(2));

    Nurbs2D_plate.cPoints{p1} = reshape( ...
        Nurbs2D_plate.nodesDeformed{p1}.', ...
        4,Nurbs2D_plate.number{p1}(1),Nurbs2D_plate.number{p1}(2));

    Nurbs2D_plate.cPoints{p2} = reshape( ...
        Nurbs2D_plate.nodesDeformed{p2}.', ...
        4,Nurbs2D_plate.number{p2}(1),Nurbs2D_plate.number{p2}(2));
end
end
