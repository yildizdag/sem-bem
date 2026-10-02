%% IP-1: CURVED NURBS SURFACE -> DELTA-CONTROLLED SEM MESH
%
% TÜBİTAK 3501 - IP-1
%
% This example creates a genuinely curved NURBS surface directly from
% control points, weights and knot vectors.
%
% The NURBS geometry and SEM discretization are independent:
%
%   1. The exact geometry is defined by the NURBS surface.
%   2. NURBS knot spans define the coarse elements.
%   3. The SEM polynomial order is selected from the delta criterion.
%   4. Chebyshev-Gauss-Lobatto points are evaluated directly on the
%      exact NURBS surface.
%
% IMPORTANT:
%
% Even for a low SEM polynomial order, the displayed element boundaries
% remain exact curved NURBS curves.
%
% Delta selection follows the previous IP-1 implementation:
%
%   reference length = global x extent
%   edge length      = straight distance between element corners
%
%   delta = max Lobatto spacing / sqrt(referenceLength*edgeLength)
%
% samePolynomialDegree = true:
%
%   The same polynomial order is enforced in xi and eta.
%   Since the mesh is connected and conforming, the largest required
%   order propagates through the entire mesh.
%
% -------------------------------------------------------------------------

clear;
clc;
close all;

%% ========================================================================
%  USER SETTINGS
%  ========================================================================

% Maximum admissible delta
delta0 = 0.15;

% Minimum SEM polynomial degree
%
% Number of Lobatto points:
%
%       N = pSEM + 1
%
minPolynomialOrder = 2;

% true:
%       Same polynomial degree in xi and eta.
%
% false:
%       Xi and eta may have different polynomial degrees.
%
samePolynomialDegree = true;

% Plot resolution for exact NURBS geometry
%
% This has NO influence on the SEM discretization.
%
nurbsPlotResolution = 61;

% Merge tolerance for coincident SEM nodes
mergeTolerance = 1e-10;

%% ========================================================================
%  CREATE CURVED NURBS SURFACE
%  ========================================================================

nurbs = create_curved_nurbs_surface();

%% ========================================================================
%  BUILD DELTA-CONTROLLED SEM MESH
%  ========================================================================

mesh = build_sem_mesh( ...
    nurbs, ...
    delta0, ...
    minPolynomialOrder, ...
    samePolynomialDegree, ...
    mergeTolerance);

%% ========================================================================
%  OUTPUT
%  ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CURVED NURBS -> DELTA-CONTROLLED SEM MESH\n');
fprintf('============================================================\n\n');

fprintf('NURBS degree in u = %d\n',nurbs.p);
fprintf('NURBS degree in v = %d\n',nurbs.q);

fprintf('Number of control points = %d x %d\n', ...
    nurbs.nu,nurbs.nv);

fprintf('Number of NURBS elements = %d\n',mesh.nel);

fprintf('\ndelta0 = %.6g\n',delta0);

fprintf('Reference x length = %.6g\n', ...
    mesh.referenceLength);

fprintf('Minimum SEM polynomial order = %d\n', ...
    minPolynomialOrder);

if samePolynomialDegree

    fprintf('Polynomial-degree mode: SAME in xi and eta\n');

    fprintf('Final common element order = %d\n', ...
        mesh.commonElementOrder);

else

    fprintf('Polynomial-degree mode: FREE\n');

    fprintf('Minimum element order = %d\n', ...
        min(mesh.polynomialOrder(:)));

    fprintf('Maximum element order = %d\n', ...
        max(mesh.polynomialOrder(:)));

end

fprintf('Number of unique SEM nodes = %d\n', ...
    size(mesh.nodes,1));

fprintf('Maximum selected delta = %.6g\n', ...
    max(mesh.edgeDelta));

fprintf('\n');

%% ========================================================================
%  ELEMENT INFORMATION
%  ========================================================================

Element = (1:mesh.nel).';

Nxi = mesh.polynums(:,1);
Neta = mesh.polynums(:,2);

DegreeXi = Nxi-1;
DegreeEta = Neta-1;

DeltaXi = mesh.elementDelta(:,1);
DeltaEta = mesh.elementDelta(:,2);

elementTable = table( ...
    Element, ...
    Nxi, ...
    Neta, ...
    DegreeXi, ...
    DegreeEta, ...
    DeltaXi, ...
    DeltaEta);

disp(elementTable);

%% ========================================================================
%  FIGURE 1
%
%  EXACT NURBS SURFACE + CONTROL NET
%  ========================================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 900 700]);

ax1 = axes;
hold(ax1,'on');

%% Exact NURBS surface

[X,Y,Z] = evaluate_complete_surface( ...
    nurbs, ...
    nurbsPlotResolution);

surf( ...
    ax1, ...
    X,Y,Z, ...
    'FaceColor',[0.82 0.90 1.00], ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0);

%% Control net

for j = 1:nurbs.nv

    P = squeeze(nurbs.controlPoints(:, :, j)).';

    plot3( ...
        ax1, ...
        P(:,1), ...
        P(:,2), ...
        P(:,3), ...
        '--', ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',1.0);

end

for i = 1:nurbs.nu

    P = squeeze(nurbs.controlPoints(:, i, :)).';

    plot3( ...
        ax1, ...
        P(:,1), ...
        P(:,2), ...
        P(:,3), ...
        '--', ...
        'Color',[0.35 0.35 0.35], ...
        'LineWidth',1.0);

end

%% Control points

CP = reshape( ...
    permute(nurbs.controlPoints,[2 3 1]), ...
    [],3);

plot3( ...
    ax1, ...
    CP(:,1), ...
    CP(:,2), ...
    CP(:,3), ...
    'o', ...
    'MarkerSize',6, ...
    'MarkerFaceColor',[0.85 0.20 0.20], ...
    'MarkerEdgeColor','k');

%% Plot exact NURBS element boundaries

plot_nurbs_element_boundaries( ...
    ax1, ...
    nurbs, ...
    mesh, ...
    nurbsPlotResolution, ...
    [0 0 0], ...
    1.0);

axis(ax1,'equal');
axis(ax1,'off');

view(ax1,[-38 27]);

title( ...
    ax1, ...
    'Curved NURBS Surface and Control Net', ...
    'FontWeight','normal');

%% ========================================================================
%  FIGURE 2
%
%  EXACT NURBS SURFACE + SEM POINTS
%  ========================================================================

figure( ...
    'Color','w', ...
    'Position',[150 100 900 700]);

ax2 = axes;
hold(ax2,'on');

%% Exact NURBS surface

surf( ...
    ax2, ...
    X,Y,Z, ...
    'FaceColor',[0.82 1.00 0.82], ...
    'EdgeColor','none', ...
    'FaceAlpha',1.0);

%% Exact NURBS element boundaries

plot_nurbs_element_boundaries( ...
    ax2, ...
    nurbs, ...
    mesh, ...
    nurbsPlotResolution, ...
    [0 0 0], ...
    1.3);

%% SEM points

plot3( ...
    ax2, ...
    mesh.nodes(:,1), ...
    mesh.nodes(:,2), ...
    mesh.nodes(:,3), ...
    '.', ...
    'Color',[0 0.20 0.85], ...
    'MarkerSize',16);

axis(ax2,'equal');
axis(ax2,'off');

view(ax2,[-38 27]);

if samePolynomialDegree

    title( ...
        ax2, ...
        sprintf( ...
        'SEM-BEM Points: \\delta_0 = %.2f, Element Order = %d', ...
        delta0, ...
        mesh.commonElementOrder), ...
        'FontWeight','normal');

else

    minOrder = min(mesh.polynomialOrder(:));
    maxOrder = max(mesh.polynomialOrder(:));

    title( ...
        ax2, ...
        sprintf( ...
        'SEM-BEM Points: \\delta_0 = %.2f, Element Orders = %d-%d', ...
        delta0, ...
        minOrder, ...
        maxOrder), ...
        'FontWeight','normal');

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  CREATE A DOUBLY-CURVED NURBS SURFACE
%  ========================================================================

function nurbs = create_curved_nurbs_surface()

    %% --------------------------------------------------------------------
    %  NURBS degree
    %
    %  Cubic in both directions
    %  --------------------------------------------------------------------

    p = 3;
    q = 3;

    %% --------------------------------------------------------------------
    %  Knot vectors
    %
    %  Four knot spans in each direction:
    %
    %       [0,0.25]
    %       [0.25,0.50]
    %       [0.50,0.75]
    %       [0.75,1.00]
    %
    %  Therefore:
    %
    %       4 x 4 = 16 NURBS elements
    %
    %  --------------------------------------------------------------------

    U = [ ...
        0 0 0 0 ...
        0.25 0.50 0.75 ...
        1 1 1 1];

    V = [ ...
        0 0 0 0 ...
        0.25 0.50 0.75 ...
        1 1 1 1];

    %% --------------------------------------------------------------------
    %  Number of control points
    %  --------------------------------------------------------------------

    nu = length(U)-p-1;
    nv = length(V)-q-1;

    assert(nu == 7);
    assert(nv == 7);

    %% --------------------------------------------------------------------
    %  Geometric dimensions
    %
    %  The surface is intentionally shallow so that it is suitable as
    %  a representative curved shell/plate geometry.
    %  --------------------------------------------------------------------

    Lx = 1.00;
    Ly = 1.00;

    % Maximum control-net rise
    H = 0.18;

    %% --------------------------------------------------------------------
    %  Control-point coordinates
    %  --------------------------------------------------------------------

    xCP = linspace(0,Lx,nu);
    yCP = linspace(0,Ly,nv);

    controlPoints = zeros(3,nu,nv);

    %% --------------------------------------------------------------------
    %  Construct a doubly-curved control net
    %
    %  z is defined ONLY at the control points.
    %
    %  The physical surface itself is subsequently generated entirely
    %  by the NURBS basis functions. Thus, z=f(x,y) is NOT used to place
    %  the SEM points.
    %
    %  A shallow dome-like geometry is used.
    %  --------------------------------------------------------------------

    for i = 1:nu

        for j = 1:nv

            x = xCP(i);
            y = yCP(j);

            xn = 2*x/Lx-1;
            yn = 2*y/Ly-1;

            % Smooth doubly-curved control net
            %
            % Boundary is lower and the central region is raised.
            %
            z = H*(1-0.55*xn^2-0.55*yn^2);

            % Additional mild coupling gives a genuinely non-cylindrical
            % doubly-curved surface.
            z = z + 0.025*xn*yn;

            controlPoints(:,i,j) = [x;y;z];

        end

    end

    %% --------------------------------------------------------------------
    %  NURBS weights
    %
    %  Non-uniform positive weights are used so that the geometry is
    %  genuinely rational rather than simply a B-spline surface.
    %
    %  The variation is intentionally mild.
    %  --------------------------------------------------------------------

    weights = ones(nu,nv);

    for i = 1:nu

        for j = 1:nv

            xn = 2*(i-1)/(nu-1)-1;
            yn = 2*(j-1)/(nv-1)-1;

            weights(i,j) = ...
                1.0 + 0.12*(1-xn^2)*(1-yn^2);

        end

    end

    %% --------------------------------------------------------------------
    %  Store
    %  --------------------------------------------------------------------

    nurbs.p = p;
    nurbs.q = q;

    nurbs.U = U;
    nurbs.V = V;

    nurbs.nu = nu;
    nurbs.nv = nv;

    nurbs.controlPoints = controlPoints;
    nurbs.weights = weights;

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  BUILD DELTA-CONTROLLED SEM MESH
%  ========================================================================

function mesh = build_sem_mesh( ...
        nurbs,delta0,minDegree,sameDegree,tol)

    %% --------------------------------------------------------------------
    %  Non-zero knot spans
    %  --------------------------------------------------------------------

    uBreaks = unique(nurbs.U);
    vBreaks = unique(nurbs.V);

    uBreaks = uBreaks(diff([uBreaks inf])>0);
    vBreaks = vBreaks(diff([vBreaks inf])>0);

    % Unique already removes repeated knots.
    %
    % Number of non-zero spans:
    %
    nElemU = numel(uBreaks)-1;
    nElemV = numel(vBreaks)-1;

    ne = nElemU*nElemV;

    %% --------------------------------------------------------------------
    %  Element definition
    %  --------------------------------------------------------------------

    elements = repmat( ...
        struct( ...
        'u',[], ...
        'v',[], ...
        'iu',0, ...
        'iv',0), ...
        ne,1);

    corners = zeros(4*ne,3);

    e = 0;

    for ev = 1:nElemV

        for eu = 1:nElemU

            e = e+1;

            ua = uBreaks(eu);
            ub = uBreaks(eu+1);

            va = vBreaks(ev);
            vb = vBreaks(ev+1);

            %% Knot-span index

            iu = find( ...
                nurbs.U <= ua, ...
                1,'last');

            if ua == nurbs.U(end)
                iu = length(nurbs.U)-nurbs.p-1;
            end

            while iu < length(nurbs.U) ...
                    && nurbs.U(iu+1) <= ua
                iu = iu+1;
            end

            iv = find( ...
                nurbs.V <= va, ...
                1,'last');

            if va == nurbs.V(end)
                iv = length(nurbs.V)-nurbs.q-1;
            end

            while iv < length(nurbs.V) ...
                    && nurbs.V(iv+1) <= va
                iv = iv+1;
            end

            % MATLAB indexing correction for span convention
            iu = find_knot_span( ...
                nurbs.nu-1, ...
                nurbs.p, ...
                0.5*(ua+ub), ...
                nurbs.U);

            iv = find_knot_span( ...
                nurbs.nv-1, ...
                nurbs.q, ...
                0.5*(va+vb), ...
                nurbs.V);

            elements(e).u = [ua ub];
            elements(e).v = [va vb];

            elements(e).iu = iu;
            elements(e).iv = iv;

            %% ------------------------------------------------------------
            %  Physical element corners
            %  ------------------------------------------------------------

            uv = [ ...
                ua va;
                ub va;
                ub vb;
                ua vb];

            for a = 1:4

                P = nurbs_surface_point( ...
                    nurbs, ...
                    uv(a,1), ...
                    uv(a,2));

                corners(4*(e-1)+a,:) = P.';

            end

        end

    end

    %% --------------------------------------------------------------------
    %  Reference x length
    %  --------------------------------------------------------------------

    sampleN = 51;

    uu = linspace( ...
        nurbs.U(nurbs.p+1), ...
        nurbs.U(end-nurbs.p), ...
        sampleN);

    vv = linspace( ...
        nurbs.V(nurbs.q+1), ...
        nurbs.V(end-nurbs.q), ...
        sampleN);

    geometryPoints = zeros(sampleN*sampleN,3);

    a = 0;

    for i = 1:sampleN

        for j = 1:sampleN

            a = a+1;

            P = nurbs_surface_point( ...
                nurbs, ...
                uu(i), ...
                vv(j));

            geometryPoints(a,:) = P.';

        end

    end

    Lref = ...
        max(geometryPoints(:,1)) ...
        - min(geometryPoints(:,1));

    assert(Lref > 0);

    %% --------------------------------------------------------------------
    %  Identify common corner nodes
    %  --------------------------------------------------------------------

    [cornerNodes,~,cornerIds] = ...
        uniquetol( ...
        corners, ...
        tol, ...
        'ByRows',true);

    elementCorners = ...
        reshape(cornerIds,4,ne).';

    %% --------------------------------------------------------------------
    %  Element edges
    %
    %       4 -------- 3
    %       |          |
    %       |          |
    %       1 -------- 2
    %
    %  --------------------------------------------------------------------

    edgeOrder = [ ...
        1 2;
        2 3;
        4 3;
        1 4];

    rawEdges = zeros(4*ne,2);

    for a = 1:4

        rawEdges((a-1)*ne+(1:ne),:) = ...
            elementCorners(:,edgeOrder(a,:));

    end

    %% --------------------------------------------------------------------
    %  Global edge IDs
    %  --------------------------------------------------------------------

    [edges,~,edgeIds] = ...
        unique( ...
        sort(rawEdges,2), ...
        'rows');

    elementEdges = ...
        reshape(edgeIds,ne,4);

    %% --------------------------------------------------------------------
    %  Polynomial-order coupling
    %  --------------------------------------------------------------------

    if sameDegree

        % Couple all edges of every element.
        %
        % Shared global edges then propagate the same degree through
        % the complete connected mesh.

        pairs = [ ...
            elementEdges(:,[1 2]);
            elementEdges(:,[1 3]);
            elementEdges(:,[1 4])];

    else

        % Couple only opposite edges.
        %
        % xi and eta may therefore have different orders.

        pairs = [ ...
            elementEdges(:,[1 3]);
            elementEdges(:,[2 4])];

    end

    G = graph( ...
        pairs(:,1), ...
        pairs(:,2), ...
        [], ...
        size(edges,1));

    groups = conncomp(G).';

    %% --------------------------------------------------------------------
    %  Physical chord length of every coarse edge
    %  --------------------------------------------------------------------

    edgeLength = vecnorm( ...
        cornerNodes(edges(:,2),:) ...
        - cornerNodes(edges(:,1),:), ...
        2,2);

    assert( ...
        all(edgeLength > 0), ...
        'Degenerate coarse edge.');

    %% --------------------------------------------------------------------
    %  Longest edge in each connected edge group
    %  --------------------------------------------------------------------

    groupLength = ...
        accumarray( ...
        groups, ...
        edgeLength, ...
        [], ...
        @max);

    groupN = zeros(size(groupLength));
    groupDelta = zeros(size(groupLength));

    %% --------------------------------------------------------------------
    %  Determine SEM order from delta
    %  --------------------------------------------------------------------

    for g = 1:numel(groupLength)

        N = minDegree+1;

        d = ip1_delta( ...
            groupLength(g), ...
            Lref, ...
            N);

        while d > delta0

            N = N+1;

            assert( ...
                N <= 10000, ...
                'delta0 requires more than 10000 points.');

            d = ip1_delta( ...
                groupLength(g), ...
                Lref, ...
                N);

        end

        groupN(g) = N;
        groupDelta(g) = d;

    end

    %% --------------------------------------------------------------------
    %  Edge number of Lobatto points
    %  --------------------------------------------------------------------

    edgeN = groupN(groups);

    %% --------------------------------------------------------------------
    %  Element Nxi and Neta
    %  --------------------------------------------------------------------

    polynums = [ ...
        edgeN(elementEdges(:,1)), ...
        edgeN(elementEdges(:,2))];

    %% --------------------------------------------------------------------
    %  Same-degree verification
    %  --------------------------------------------------------------------

    if sameDegree

        assert( ...
            all(polynums(:,1) == polynums(:,2)), ...
            'Nxi and Neta must be identical.');

        commonOrders = ...
            unique(polynums(:)-1);

        assert( ...
            numel(commonOrders) == 1, ...
            'Connected mesh contains multiple orders.');

        commonElementOrder = ...
            commonOrders(1);

    else

        commonElementOrder = [];

    end

    %% --------------------------------------------------------------------
    %  Generate SEM points on exact NURBS geometry
    %  --------------------------------------------------------------------

    elementNodes = cell(ne,1);

    allNodes = [];

    for e = 1:ne

        nx = polynums(e,1);
        ny = polynums(e,2);

        xi = lobatto_points(nx);
        eta = lobatto_points(ny);

        elementXYZ = zeros(nx*ny,3);

        a = 0;

        for i = 1:nx

            for j = 1:ny

                a = a+1;

                %% Map Lobatto coordinate [-1,1]
                %  to NURBS knot span

                u = ...
                    (1-xi(i))*elements(e).u(1)/2 ...
                    +(1+xi(i))*elements(e).u(2)/2;

                v = ...
                    (1-eta(j))*elements(e).v(1)/2 ...
                    +(1+eta(j))*elements(e).v(2)/2;

                %% Exact NURBS evaluation

                P = nurbs_surface_point( ...
                    nurbs,u,v);

                elementXYZ(a,:) = P.';

            end

        end

        startID = size(allNodes,1);

        allNodes = [allNodes;elementXYZ]; %#ok<AGROW>

        elementNodes{e} = ...
            startID+(1:size(elementXYZ,1));

    end

    %% --------------------------------------------------------------------
    %  Merge coincident SEM nodes
    %  --------------------------------------------------------------------

    [nodes,~,globalIds] = ...
        uniquetol( ...
        allNodes, ...
        tol, ...
        'ByRows',true);

    nodeConn = cell(ne,1);

    for e = 1:ne

        nodeConn{e} = ...
            globalIds(elementNodes{e}).';

    end

    %% --------------------------------------------------------------------
    %  Element delta
    %  --------------------------------------------------------------------

    edgeDelta = ...
        groupDelta(groups);

    elementDelta = [ ...
        edgeDelta(elementEdges(:,1)), ...
        edgeDelta(elementEdges(:,2))];

    %% --------------------------------------------------------------------
    %  Store
    %  --------------------------------------------------------------------

    mesh.nel = ne;

    mesh.elements = elements;

    mesh.referenceLength = Lref;

    mesh.cornerNodes = cornerNodes;
    mesh.elementCorners = elementCorners;

    mesh.edges = edges;
    mesh.elementEdges = elementEdges;

    mesh.edgesgroup = groups;

    mesh.edgeLength = edgeLength;
    mesh.edgeDelta = edgeDelta;

    mesh.polynums = polynums;
    mesh.polynomialOrder = polynums-1;

    mesh.samePolynomialDegree = sameDegree;

    mesh.commonElementOrder = ...
        commonElementOrder;

    mesh.nodes = nodes;
    mesh.nodeConn = nodeConn;

    mesh.elementDelta = elementDelta;

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  DELTA EVALUATION
%
%  Same definition as the previous IP-1 code.
%  ========================================================================

function d = ip1_delta(edgeLength,referenceLength,N)

    x = ...
        (edgeLength/2)*lobatto_points(N);

    d = ...
        max(diff(x)) ...
        /sqrt(referenceLength*edgeLength);

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  CHEBYSHEV-GAUSS-LOBATTO POINTS
%
%  Returns ascending coordinates from -1 to 1.
%  ========================================================================

function x = lobatto_points(N)

    assert(N >= 2);

    k = (0:N-1).';

    x = ...
        -cos(pi*k/(N-1));

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  EVALUATE COMPLETE EXACT NURBS SURFACE
%  ========================================================================

function [X,Y,Z] = evaluate_complete_surface(nurbs,nPlot)

    u = linspace( ...
        nurbs.U(nurbs.p+1), ...
        nurbs.U(end-nurbs.p), ...
        nPlot);

    v = linspace( ...
        nurbs.V(nurbs.q+1), ...
        nurbs.V(end-nurbs.q), ...
        nPlot);

    X = zeros(nPlot,nPlot);
    Y = zeros(nPlot,nPlot);
    Z = zeros(nPlot,nPlot);

    for i = 1:nPlot

        for j = 1:nPlot

            P = nurbs_surface_point( ...
                nurbs,u(i),v(j));

            X(j,i) = P(1);
            Y(j,i) = P(2);
            Z(j,i) = P(3);

        end

    end

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  PLOT EXACT NURBS ELEMENT BOUNDARIES
%  ========================================================================

function plot_nurbs_element_boundaries( ...
        ax,nurbs,mesh,nPlot,lineColor,lineWidth)

    for e = 1:mesh.nel

        ua = mesh.elements(e).u(1);
        ub = mesh.elements(e).u(2);

        va = mesh.elements(e).v(1);
        vb = mesh.elements(e).v(2);

        t = linspace(0,1,nPlot);

        boundary = zeros(4*nPlot,3);

        %% Bottom edge

        for a = 1:nPlot

            u = ua+t(a)*(ub-ua);
            v = va;

            P = nurbs_surface_point(nurbs,u,v);

            boundary(a,:) = P.';

        end

        %% Right edge

        for a = 1:nPlot

            u = ub;
            v = va+t(a)*(vb-va);

            P = nurbs_surface_point(nurbs,u,v);

            boundary(nPlot+a,:) = P.';

        end

        %% Top edge

        for a = 1:nPlot

            u = ub-t(a)*(ub-ua);
            v = vb;

            P = nurbs_surface_point(nurbs,u,v);

            boundary(2*nPlot+a,:) = P.';

        end

        %% Left edge

        for a = 1:nPlot

            u = ua;
            v = vb-t(a)*(vb-va);

            P = nurbs_surface_point(nurbs,u,v);

            boundary(3*nPlot+a,:) = P.';

        end

        plot3( ...
            ax, ...
            boundary(:,1), ...
            boundary(:,2), ...
            boundary(:,3), ...
            '-', ...
            'Color',lineColor, ...
            'LineWidth',lineWidth);

    end

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  EXACT NURBS SURFACE POINT
%  ========================================================================

function P = nurbs_surface_point(nurbs,u,v)

    p = nurbs.p;
    q = nurbs.q;

    U = nurbs.U;
    V = nurbs.V;

    n = nurbs.nu-1;
    m = nurbs.nv-1;

    %% Knot spans

    spanU = find_knot_span(n,p,u,U);
    spanV = find_knot_span(m,q,v,V);

    %% B-spline basis functions

    Nu = basis_functions(spanU,u,p,U);
    Nv = basis_functions(spanV,v,q,V);

    %% Rational surface

    numerator = zeros(3,1);
    denominator = 0;

    for a = 0:p

        i = spanU-p+a;

        for b = 0:q

            j = spanV-q+b;

            w = nurbs.weights(i+1,j+1);

            B = ...
                Nu(a+1) ...
                *Nv(b+1) ...
                *w;

            numerator = ...
                numerator ...
                + B*nurbs.controlPoints(:,i+1,j+1);

            denominator = ...
                denominator+B;

        end

    end

    assert( ...
        abs(denominator) > eps, ...
        'Zero NURBS denominator.');

    P = numerator/denominator;

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  FIND KNOT SPAN
%
%  Piegl & Tiller convention converted to MATLAB implementation.
%  Returned span is ZERO-BASED.
%  ========================================================================

function span = find_knot_span(n,p,u,U)

    %% Right boundary

    if u >= U(n+2)

        span = n;
        return;

    end

    %% Left boundary

    if u <= U(p+1)

        span = p;
        return;

    end

    low = p;
    high = n+1;

    mid = floor((low+high)/2);

    while ...
            u < U(mid+1) ...
            || u >= U(mid+2)

        if u < U(mid+1)

            high = mid;

        else

            low = mid;

        end

        mid = floor((low+high)/2);

    end

    span = mid;

end

%% ========================================================================
%  LOCAL FUNCTION
%
%  NON-ZERO B-SPLINE BASIS FUNCTIONS
%
%  Piegl & Tiller Algorithm A2.2
%
%  span is ZERO-BASED.
%  ========================================================================

function N = basis_functions(span,u,p,U)

    N = zeros(1,p+1);

    left = zeros(1,p+1);
    right = zeros(1,p+1);

    N(1) = 1.0;

    for j = 1:p

        left(j+1) = ...
            u-U(span+2-j);

        right(j+1) = ...
            U(span+j+1)-u;

        saved = 0.0;

        for r = 0:j-1

            denominator = ...
                right(r+2) ...
                + left(j-r+1);

            if abs(denominator) < eps

                temp = 0;

            else

                temp = ...
                    N(r+1)/denominator;

            end

            N(r+1) = ...
                saved ...
                + right(r+2)*temp;

            saved = ...
                left(j-r+1)*temp;

        end

        N(j+1) = saved;

    end

end