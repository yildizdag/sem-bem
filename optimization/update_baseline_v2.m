function [Nurbs2D_plate,Nurbs2D_stiff] = update_baseline_v2(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp)
%UPDATE_BASELINE_V2 Updates plate-stiffener geometry by moving CP groups.
%
% For every design variable, one complete control-point group is moved:
%   1) interface CP on plate patch 1,
%   2) adjacent inward CP on plate patch 1,
%   3) coincident interface CP on plate patch 2,
%   4) adjacent inward CP on plate patch 2,
%   5) all stiffener CPs lying on the same physical line as the interface CP.
%
% The correspondence between the two plate edges is determined from the
% physical coordinates. Therefore, opposite edge orientations do not cause
% incorrect pairing.
%
% Current design convention:
%   dcp(g) = displacement of group g in the global x direction.
%
% The original geometry is used at every call, so geometry changes do not
% accumulate between objective-function evaluations.

tol = 1.0e-8;

%% Copy baseline control points
Nurbs2D_plate.nodesDeformed = Nurbs2D_plate.nodes;
Nurbs2D_stiff.nodesDeformed = Nurbs2D_stiff.nodes;

%% Number of design variables
nvars = 0;
for c = 1:size(pconn,1)
    p1 = pconn(c,1);
    e1 = pconn(c,2);
    nvars = nvars + size(Nurbs2D_plate.movingCP{p1,e1},2);
end

dcp = dcp(:).';

if numel(dcp) ~= nvars
    error('update_baseline_v2:InvalidDesignVector', ...
        'dcp must contain exactly %d values.',nvars);
end

%% Move every plate-stiffener group
offset = 0;

for c = 1:size(pconn,1)

    p1 = pconn(c,1);
    e1 = pconn(c,2);

    p2 = pconn(c,3);
    e2 = pconn(c,4);

    s  = pconn(c,5);

    moving1 = Nurbs2D_plate.movingCP{p1,e1};
    moving2 = Nurbs2D_plate.movingCP{p2,e2};

    ngroup = size(moving1,2);

    if size(moving1,1) < 2 || size(moving2,1) < 2
        error('update_baseline_v2:MissingAdjacentCP', ...
            'Each plate movingCP array must contain interface and adjacent rows.');
    end

    for g = 1:ngroup

        dx = dcp(offset+g);
        dy = 0.0;
        dz = 0.0;

        % -------------------------------------------------------------
        % Plate patch 1: interface CP and adjacent inward CP
        % -------------------------------------------------------------
        cp1_interface = moving1(1,g);
        cp1_adjacent  = moving1(2,g);

        xyzInterface = Nurbs2D_plate.nodes{p1}(cp1_interface,1:3);

        % -------------------------------------------------------------
        % Plate patch 2: locate the physically coincident interface CP.
        % Do not assume that both connected edges have the same direction.
        % -------------------------------------------------------------
        interfaceCandidates2 = moving2(1,:);

        delta2 = Nurbs2D_plate.nodes{p2}(interfaceCandidates2,1:3) ...
               - xyzInterface;

        [minDist2,g2] = min(vecnorm(delta2,2,2));

        if minDist2 > tol
            error('update_baseline_v2:PlateInterfaceMismatch', ...
                ['No coincident control point was found on plate patch %d ' ...
                 'for interface CP %d of plate patch %d.'], ...
                 p2,cp1_interface,p1);
        end

        cp2_interface = moving2(1,g2);
        cp2_adjacent  = moving2(2,g2);

        % -------------------------------------------------------------
        % Stiffener: all CPs on the same physical line.
        % For a vertical stiffener, these CPs have the same x-y position
        % and may have different z coordinates.
        % -------------------------------------------------------------
        stiffGroup = [];

        if s > 0
            deltaXY = Nurbs2D_stiff.nodes{s}(:,1:2) ...
                    - xyzInterface(1:2);

            stiffGroup = find(vecnorm(deltaXY,2,2) <= tol);

            if isempty(stiffGroup)
                error('update_baseline_v2:StiffenerLineNotFound', ...
                    ['No stiffener control points were found on the line ' ...
                     'passing through plate interface CP %d of patch %d.'], ...
                     cp1_interface,p1);
            end
        end

        % -------------------------------------------------------------
        % Move the complete group in exactly the same direction
        % -------------------------------------------------------------
        groupDisp = [dx,dy,dz];

        plateGroup1 = unique([cp1_interface;cp1_adjacent]);
        plateGroup2 = unique([cp2_interface;cp2_adjacent]);

        Nurbs2D_plate.nodesDeformed{p1}(plateGroup1,1:3) = ...
            Nurbs2D_plate.nodes{p1}(plateGroup1,1:3) + groupDisp;

        Nurbs2D_plate.nodesDeformed{p2}(plateGroup2,1:3) = ...
            Nurbs2D_plate.nodes{p2}(plateGroup2,1:3) + groupDisp;

        if ~isempty(stiffGroup)
            Nurbs2D_stiff.nodesDeformed{s}(stiffGroup,1:3) = ...
                Nurbs2D_stiff.nodes{s}(stiffGroup,1:3) + groupDisp;
        end
    end

    offset = offset + ngroup;
end

%% Rebuild plate NURBS control nets
for p = 1:Nurbs2D_plate.numpatch
    Nu = Nurbs2D_plate.number{p}(1);
    Nv = Nurbs2D_plate.number{p}(2);

    Nurbs2D_plate.cPoints{p} = reshape( ...
        Nurbs2D_plate.nodesDeformed{p}.',4,Nu,Nv);
end

%% Rebuild stiffener NURBS control nets
for s = 1:Nurbs2D_stiff.numpatch
    Nu = Nurbs2D_stiff.number{s}(1);
    Nv = Nurbs2D_stiff.number{s}(2);

    Nurbs2D_stiff.cPoints{s} = reshape( ...
        Nurbs2D_stiff.nodesDeformed{s}.',4,Nu,Nv);
end

end
