%% IP-1: NURBS geometry -> delta-controlled SEM mesh
%
% Geometry evaluation follows sem_core/sem2Dmesh.m.
%
% Delta selection:
%   - reference length = global x extent;
%   - edge length = straight distance between its corner points;
%   - shared edges have a common sample count;
%   - each edge group uses its longest edge;
%   - increase the sample count by one until delta <= delta0.
%
% Polynomial-degree mode:
%
%   samePolynomialDegree = true
%       Enforces the same SEM polynomial degree in xi and eta.
%       For a connected conforming mesh, the same polynomial degree
%       propagates through all elements.
%
%   samePolynomialDegree = false
%       Xi and eta are allowed to use different polynomial degrees.
%
% IMPORTANT:
% The NURBS geometry and SEM discretization are independent.
% The exact NURBS surface is evaluated separately from the SEM points.
% Therefore, even for a low-order SEM discretization, curved NURBS
% boundaries remain curved.
%
% This example generates the SEM mesh only.
% It does not call the SEM/BEM solver.

clear;
clc;


%% ========================================================================
%  REPOSITORY AND FUNCTION PATHS
%  ========================================================================

% Folder containing this script
scriptFolder = fileparts(mfilename('fullpath'));

% Repository root
% This assumes that this script is located one folder below the repo root.
repoRoot = fileparts(scriptFolder);

% Main SEM/NURBS folder
semCoreFolder = fullfile(repoRoot,'sem_core');

assert(isfolder(semCoreFolder), ...
    'SEM core folder not found: %s',semCoreFolder);

% Add sem_core and all subfolders
addpath(genpath(semCoreFolder),'-begin');


%% ========================================================================
%  GEOMETRY PATH
%  ========================================================================

geometryFolder = fullfile(semCoreFolder,'geometry');

assert(isfolder(geometryFolder), ...
    'Geometry folder not found: %s',geometryFolder);


%% ========================================================================
%  CHECK REQUIRED FUNCTIONS
%  ========================================================================

requiredFunctions = { ...
    'iga2Dmesh', ...
    'dersbasisfuns', ...
    'derRat2DBasisFuns', ...
    'lobat', ...
    'cheb', ...
    'derivative', ...
    'InnerProduct'};

for iFunction = 1:numel(requiredFunctions)

    functionName = requiredFunctions{iFunction};

    assert(exist(functionName,'file') ~= 0, ...
        'Required function "%s" was not found on the MATLAB path.', ...
        functionName);

end


%% ========================================================================
%  USER SETTINGS
%  ========================================================================

% Maximum admissible delta
delta0 = 0.15;

% Minimum SEM polynomial degree
%
% Number of Lobatto points:
%
%       N = p + 1
%
minPolynomialOrder = 2;

% Polynomial-degree constraint
%
% true:
%       Same polynomial order in xi and eta.
%       For a connected conforming mesh, all elements consequently
%       obtain the same polynomial order.
%
% false:
%       DegreeXi and DegreeEta may differ.
%
samePolynomialDegree = true;

% Geometry file prefix
fileName = 'plate_cutout_';

% Number of NURBS patches
numPatch = 4;

% Number of shell DOFs per node
shell_dof = 6;

% Relative tolerance used when merging coincident SEM nodes
mergeTolerance = 1e-5;

% Plot figures
showFigures = true;

% Dense sampling used ONLY to display the exact NURBS geometry.
% This parameter has no influence on the SEM discretization.
nurbsPlotResolution = 41;


%% ========================================================================
%  GEOMETRY FILE PREFIX
%  ========================================================================

geometryPrefix = fullfile(geometryFolder,fileName);


%% ========================================================================
%  BUILD SEM MESH
%  ========================================================================

[Nurbs2D,sem2D,elementSummary,edgeSummary] = ip1_build_mesh( ...
    geometryPrefix, ...
    numPatch, ...
    delta0, ...
    minPolynomialOrder, ...
    samePolynomialDegree, ...
    shell_dof, ...
    mergeTolerance);


%% ========================================================================
%  OUTPUT
%  ========================================================================

fprintf('\n');

fprintf('Geometry: %s (%d patches)\n', ...
    fileName,numPatch);

fprintf('delta0 = %.6g\n', ...
    delta0);

fprintf('Minimum polynomial degree = %d (%d Lobatto points)\n', ...
    minPolynomialOrder,minPolynomialOrder+1);

if samePolynomialDegree

    fprintf( ...
        'Polynomial-degree mode: SAME in xi and eta.\n');

    % All elements should have the same polynomial order
    elementOrders = unique(sem2D.polynomialOrder(:));

    assert(numel(elementOrders) == 1, ...
        ['samePolynomialDegree=true, but more than one element ' ...
         'polynomial order was detected.']);

    fprintf('Final element polynomial order = %d\n', ...
        elementOrders);

else

    fprintf( ...
        'Polynomial-degree mode: FREE in xi and eta.\n');

end

fprintf('Reference x length = %.6g\n', ...
    sem2D.referenceLength);

fprintf( ...
    '%d elements, %d edge groups, %d unique SEM nodes, %d DOFs\n', ...
    sem2D.nel, ...
    max(sem2D.edgesgroup), ...
    sem2D.non, ...
    sem2D.dof);

fprintf( ...
    'Maximum selected delta = %.6g\n', ...
    max(sem2D.edgeDelta));

fprintf( ...
    'Minimum sampled surface J = %.6g\n', ...
    min(cellfun(@min,sem2D.J)));

fprintf('\n');

disp(elementSummary);

fprintf([ ...
    '\nWorkspace variables:\n' ...
    '  Nurbs2D\n' ...
    '  sem2D\n' ...
    '  elementSummary\n' ...
    '  edgeSummary\n\n' ...
    'sem2D.polynums(e,:) = [Nxi Neta]\n' ...
    'SEM polynomial degree = polynums - 1\n\n']);


%% ========================================================================
%  PLOT
%  ========================================================================

if showFigures

    ip1_plot_mesh( ...
        Nurbs2D, ...
        sem2D, ...
        nurbsPlotResolution);

end


%% ========================================================================
%  LOCAL FUNCTION:
%  BUILD DELTA-CONTROLLED SEM MESH
%  ========================================================================

function [nurbs,mesh,summary,edgeSummary] = ip1_build_mesh( ...
        prefix,np,delta0,minDegree,sameDegree,ndof,tol)


    %% --------------------------------------------------------------------
    %  Input validation
    %  --------------------------------------------------------------------

    validateattributes( ...
        delta0, ...
        {'numeric'}, ...
        {'scalar','real','finite','positive'});

    validateattributes( ...
        minDegree, ...
        {'numeric'}, ...
        {'scalar','integer','>=',1});

    validateattributes( ...
        sameDegree, ...
        {'logical','numeric'}, ...
        {'scalar'});

    validateattributes( ...
        np, ...
        {'numeric'}, ...
        {'scalar','integer','positive'});

    validateattributes( ...
        ndof, ...
        {'numeric'}, ...
        {'scalar','integer','positive'});

    validateattributes( ...
        tol, ...
        {'numeric'}, ...
        {'scalar','real','finite','positive'});

    sameDegree = logical(sameDegree);


    %% --------------------------------------------------------------------
    %  Check geometry files
    %  --------------------------------------------------------------------

    for k = 1:np

        assert( ...
            isfile([prefix num2str(k)]), ...
            'Missing geometry file: %s%d', ...
            prefix,k);

    end


    %% --------------------------------------------------------------------
    %  Read NURBS geometry
    %  --------------------------------------------------------------------

    handlesBefore = openedFiles;

    fileCleanup = onCleanup( ...
        @() ip1_close_new_files(handlesBefore)); %#ok<NASGU>

    nurbs = iga2Dmesh(prefix,np,1);

    ip1_close_new_files(handlesBefore);


    %% --------------------------------------------------------------------
    %  Original NURBS knot-span elements
    %  --------------------------------------------------------------------

    ne = sum(cellfun(@double,nurbs.nel));

    elements = repmat( ...
        struct( ...
        'patch',0, ...
        'localElement',0, ...
        'iu',0, ...
        'iv',0, ...
        'u',[], ...
        'v',[], ...
        'CP',[], ...
        'order',[]), ...
        ne,1);

    corners = zeros(4*ne,3);

    boundary = cell(ne,1);

    geometryNodes = zeros(25*ne,3);

    e = 0;


    %% --------------------------------------------------------------------
    %  Loop over NURBS patches and knot-span elements
    %  --------------------------------------------------------------------

    for k = 1:np

        for localElement = 1:nurbs.nel{k}

            e = e+1;


            %% Element knot-span indices

            iu = nurbs.INC{k}( ...
                nurbs.IEN{k}(1,localElement),1);

            iv = nurbs.INC{k}( ...
                nurbs.IEN{k}(1,localElement),2);


            %% Parametric knot-span limits

            u = nurbs.knots.U{k}(iu:iu+1);

            v = nurbs.knots.V{k}(iv:iv+1);

            assert( ...
                u(2)>u(1) && v(2)>v(1), ...
                'Zero/negative knot span in patch %d, element %d.', ...
                k,localElement);


            %% NURBS order and control points

            order = nurbs.order{k};

            CP = nurbs.cPoints{k}( ...
                :, ...
                iu-order(1)+1:iu, ...
                iv-order(2)+1:iv);


            %% Store element information

            elements(e) = struct( ...
                'patch',k, ...
                'localElement',localElement, ...
                'iu',iu, ...
                'iv',iv, ...
                'u',u, ...
                'v',v, ...
                'CP',CP, ...
                'order',order);


            %% ------------------------------------------------------------
            %  Physical element corners
            %  ------------------------------------------------------------

            cornerParams = [ ...
                -1 -1;
                 1 -1;
                 1  1;
                -1  1];

            for a = 1:4

                s = ip1_surface( ...
                    nurbs, ...
                    elements(e), ...
                    cornerParams(a,1), ...
                    cornerParams(a,2));

                corners(4*(e-1)+a,:) = ...
                    s(:,1,1).';

            end


            %% ------------------------------------------------------------
            %  Geometry points for reference-length calculation
            %  ------------------------------------------------------------

            [x,y] = ndgrid( ...
                linspace(-1,1,5));

            for a = 1:25

                s = ip1_surface( ...
                    nurbs, ...
                    elements(e), ...
                    x(a), ...
                    y(a));

                geometryNodes(25*(e-1)+a,:) = ...
                    s(:,1,1).';

            end


            %% ------------------------------------------------------------
            %  Dense exact NURBS boundary
            %  ------------------------------------------------------------

            t = linspace(-1,1,31).';

            edgeParams = [ ...
                t,              -ones(size(t));
                ones(size(t)),   t;
                flipud(t),       ones(size(t));
                -ones(size(t)),  flipud(t)];

            boundary{e} = ...
                zeros(size(edgeParams,1),3);

            for a = 1:size(edgeParams,1)

                s = ip1_surface( ...
                    nurbs, ...
                    elements(e), ...
                    edgeParams(a,1), ...
                    edgeParams(a,2));

                boundary{e}(a,:) = ...
                    s(:,1,1).';

            end

        end

    end


    %% --------------------------------------------------------------------
    %  Reference length
    %  --------------------------------------------------------------------

    Lref = ...
        max(geometryNodes(:,1)) ...
        - min(geometryNodes(:,1));

    assert( ...
        isfinite(Lref) && Lref>0, ...
        'Delta rule requires a positive global x extent.');


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
    %  Local edges
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

    reversed = ...
        reshape( ...
        rawEdges(:,1)>rawEdges(:,2), ...
        ne,4);


    %% --------------------------------------------------------------------
    %  Edge grouping and polynomial-degree constraint
    %  --------------------------------------------------------------------

    if sameDegree

        % Couple all four edges of each element.
        %
        % Shared edges already have the same global edge ID.
        % Therefore this propagates the common polynomial order through
        % connected neighboring elements.

        pairs = [ ...
            elementEdges(:,[1 2]);
            elementEdges(:,[1 3]);
            elementEdges(:,[1 4])];

    else

        % Only opposite edges are coupled.
        %
        % Xi and eta polynomial orders can therefore differ.

        pairs = [ ...
            elementEdges(:,[1 3]);
            elementEdges(:,[2 4])];

    end


    %% --------------------------------------------------------------------
    %  Connected edge groups
    %  --------------------------------------------------------------------

    G = graph( ...
        pairs(:,1), ...
        pairs(:,2), ...
        [], ...
        size(edges,1));

    groups = conncomp(G).';


    %% --------------------------------------------------------------------
    %  Physical edge lengths
    %  --------------------------------------------------------------------

    edgeLength = vecnorm( ...
        cornerNodes(edges(:,2),:) ...
        - cornerNodes(edges(:,1),:), ...
        2,2);

    assert( ...
        all(isfinite(edgeLength) & edgeLength>0), ...
        'Degenerate coarse edge.');


    %% --------------------------------------------------------------------
    %  Longest physical edge in each edge group
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
    %  Determine SEM sampling from delta
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
                N<=10000, ...
                'delta0 requires more than 10000 points per edge.');

            d = ip1_delta( ...
                groupLength(g), ...
                Lref, ...
                N);

        end

        groupN(g) = N;

        groupDelta(g) = d;

    end


    %% --------------------------------------------------------------------
    %  Number of Lobatto points on each global edge
    %  --------------------------------------------------------------------

    edgeN = groupN(groups);


    %% --------------------------------------------------------------------
    %  Number of SEM points in xi and eta directions
    %  --------------------------------------------------------------------

    polynums = [ ...
        edgeN(elementEdges(:,1)), ...
        edgeN(elementEdges(:,2))];


    %% --------------------------------------------------------------------
    %  Verify same-degree condition
    %  --------------------------------------------------------------------

    if sameDegree

        assert( ...
            all(polynums(:,1)==polynums(:,2)), ...
            ['Internal error: same-degree mode failed ' ...
             'to produce Nxi = Neta.']);

        % Since the conforming connected mesh is intended to have one
        % common element order, verify this explicitly.

        finalOrders = unique(polynums(:,1)-1);

        assert(numel(finalOrders)==1, ...
            ['samePolynomialDegree=true, but the connected mesh ' ...
             'contains more than one polynomial order.']);

    end


    %% --------------------------------------------------------------------
    %  Initialize SEM mesh structure
    %  --------------------------------------------------------------------

    mesh.nel = ne;

    mesh.shell_dof = ndof;

    mesh.delta0 = delta0;

    mesh.minPolynomialOrder = minDegree;

    mesh.samePolynomialDegree = sameDegree;

    mesh.referenceLength = Lref;

    mesh.polynums = polynums;

    mesh.polynomialOrder = polynums-1;

    % Convenient element-order vector
    mesh.elementOrder = polynums(:,1)-1;

    if sameDegree
        mesh.commonElementOrder = mesh.elementOrder(1);
    else
        mesh.commonElementOrder = [];
    end

    mesh.cornerNodes = cornerNodes;

    mesh.elementCorners = elementCorners;

    mesh.edges = edges;

    mesh.elementedges = ...
        [elementEdges,reversed];

    mesh.edgesgroup = groups;

    mesh.edgeLength = edgeLength;

    mesh.edgeGroupLength = ...
        groupLength(groups);

    mesh.edgePolynums = edgeN;

    mesh.edgeDelta = ...
        groupDelta(groups);

    mesh.N = polynums;

    mesh.boundary = boundary;

    mesh.patch = ...
        [elements.patch].';

    mesh.localElement = ...
        [elements.localElement].';

    % Store NURBS knot-span element information for plotting
    mesh.geometryElements = elements;


    %% --------------------------------------------------------------------
    %  SEM node offsets
    %  --------------------------------------------------------------------

    offsets = [ ...
        0;
        cumsum(prod(polynums,2))];

    nodeData = ...
        zeros(offsets(end),3);


    %% --------------------------------------------------------------------
    %  Variable-size element-wise quantities
    %  --------------------------------------------------------------------

    fields = { ...
        'nodeConn', ...
        'conn', ...
        'J', ...
        'Jmat', ...
        'InvJmat', ...
        'Kappa', ...
        'R', ...
        'operators', ...
        'xi', ...
        'eta'};

    for f = 1:numel(fields)

        mesh.(fields{f}) = ...
            cell(ne,1);

    end


    %% ====================================================================
    %  EVALUATE NURBS GEOMETRY AT SEM LOBATTO POINTS
    %  ====================================================================

    for e = 1:ne

        nx = polynums(e,1);

        ny = polynums(e,2);

        nn = nx*ny;


        %% Lobatto points

        xi = lobat(nx);

        eta = lobat(ny);

        mesh.xi{e} = xi;

        mesh.eta{e} = eta;


        %% Allocate element arrays

        mesh.J{e} = ...
            zeros(nn,1);

        mesh.Jmat{e} = ...
            zeros(3,2,nn);

        mesh.InvJmat{e} = ...
            zeros(2,2,nn);

        mesh.Kappa{e} = ...
            zeros(2,nn);

        mesh.R{e} = ...
            zeros(3,3,nn);


        %% Parametric mapping factors

        du = diff(elements(e).u)/2;

        dv = diff(elements(e).v)/2;

        a = 0;


        %% ---------------------------------------------------------------
        %  Loop over SEM points
        %  ---------------------------------------------------------------

        for i = 1:nx

            for j = 1:ny

                a = a+1;


                %% Exact NURBS evaluation at Lobatto point

                s = ip1_surface( ...
                    nurbs, ...
                    elements(e), ...
                    xi(i), ...
                    eta(j));


                %% Covariant basis vectors

                A1 = s(:,2,1);

                A2 = s(:,1,2);

                A = [A1,A2];


                %% Metric tensor

                metric = A.'*A;


                %% Surface normal

                normal = cross(A1,A2);

                area = norm(normal);

                assert( ...
                    isfinite(area) ...
                    && area>0 ...
                    && rcond(metric)>eps, ...
                    'Degenerate surface: element %d, sample %d.', ...
                    e,a);

                normal = normal/area;


                %% Local orthonormal basis

                t1 = A1/norm(A1);

                t2 = cross(normal,t1);

                t2 = t2/norm(t2);


                %% Contravariant basis

                Ac = A/metric;


                %% Surface curvature

                curvature = [ ...
                    dot(s(:,3,1),normal), ...
                    dot(s(:,2,2),normal);
                    dot(s(:,2,2),normal), ...
                    dot(s(:,1,3),normal)];


                %% Physical SEM point

                nodeData(offsets(e)+a,:) = ...
                    s(:,1,1).';


                %% Surface Jacobian

                mesh.J{e}(a) = ...
                    area*du*dv;


                %% Surface Jacobian matrix

                mesh.Jmat{e}(:,:,a) = A;


                %% Inverse Jacobian

                mesh.InvJmat{e}(:,:,a) = [ ...
                    dot(t1,Ac(:,1))/du, ...
                    dot(t1,Ac(:,2))/dv;
                    dot(t2,Ac(:,1))/du, ...
                    dot(t2,Ac(:,2))/dv];


                %% Principal curvatures

                mesh.Kappa{e}(:,a) = ...
                    abs(eig(metric\curvature));


                %% Local coordinate transformation

                mesh.R{e}(:,:,a) = ...
                    [t1,t2,normal];

            end

        end


        %% SEM operators

        mesh.operators{e} = ...
            ip1_operators(nx,ny);

    end


    %% ====================================================================
    %  MERGE COINCIDENT SEM POINTS
    %  ====================================================================

    [mesh.nodes,~,ids] = ...
        uniquetol( ...
        nodeData, ...
        tol, ...
        'ByRows',true);


    %% --------------------------------------------------------------------
    %  Global connectivity
    %  --------------------------------------------------------------------

    for e = 1:ne

        nc = ids( ...
            offsets(e) ...
            +(1:prod(polynums(e,:)))).';

        mesh.nodeConn{e} = nc;

        conn = ...
            zeros(1,ndof*numel(nc));

        for d = 1:ndof

            conn(d:ndof:end) = ...
                ndof*nc-(ndof-d);

        end

        mesh.conn{e} = conn;

    end


    %% --------------------------------------------------------------------
    %  Number of global SEM nodes and DOFs
    %  --------------------------------------------------------------------

    mesh.non = ...
        size(mesh.nodes,1);

    mesh.dof = ...
        ndof*mesh.non;


    %% ====================================================================
    %  VERIFY SHARED-EDGE COMPATIBILITY
    %  ====================================================================

    edgeNodes = ...
        cell(size(edges,1),1);

    for e = 1:ne

        nx = polynums(e,1);

        ny = polynums(e,2);


        %% Local SEM edge node IDs

        localEdges = { ...
            1:ny:(nx-1)*ny+1, ...
            (nx-1)*ny+(1:ny), ...
            ny:ny:nx*ny, ...
            1:ny};


        %% Check each edge

        for a = 1:4

            current = ...
                mesh.nodeConn{e}(localEdges{a});

            if reversed(e,a)

                current = ...
                    fliplr(current);

            end

            id = ...
                elementEdges(e,a);

            if isempty(edgeNodes{id})

                edgeNodes{id} = ...
                    current;

            else

                assert( ...
                    isequal(edgeNodes{id},current), ...
                    ['Shared edge %d has incompatible sampled ' ...
                     'coordinates. Matching endpoints alone do not ' ...
                     'guarantee an identical curve/parameterization.'], ...
                    id);

            end

        end

    end

    mesh.edgeNodes = edgeNodes;


    %% ====================================================================
    %  SUMMARY TABLES
    %  ====================================================================

    deltaXi = ...
        mesh.edgeDelta(elementEdges(:,1));

    deltaEta = ...
        mesh.edgeDelta(elementEdges(:,2));


    %% Element summary

    summary = table( ...
        (1:ne).', ...
        mesh.patch, ...
        mesh.localElement, ...
        polynums(:,1), ...
        polynums(:,2), ...
        polynums(:,1)-1, ...
        polynums(:,2)-1, ...
        deltaXi, ...
        deltaEta, ...
        cellfun(@min,mesh.J), ...
        'VariableNames',{ ...
        'Element', ...
        'Patch', ...
        'PatchElement', ...
        'Nxi', ...
        'Neta', ...
        'DegreeXi', ...
        'DegreeEta', ...
        'DeltaXi', ...
        'DeltaEta', ...
        'MinSurfaceJ'});


    %% Edge summary

    edgeSummary = table( ...
        (1:size(edges,1)).', ...
        groups, ...
        edgeLength, ...
        groupLength(groups), ...
        edgeN, ...
        mesh.edgeDelta, ...
        'VariableNames',{ ...
        'Edge', ...
        'Group', ...
        'ChordLength', ...
        'GroupMaxLength', ...
        'Points', ...
        'Delta'});

end


%% ========================================================================
%  LOCAL FUNCTION:
%  EXACT NURBS SURFACE EVALUATION
%  ========================================================================

function s = ip1_surface(nurbs,el,xi,eta)

    %% Map [-1,1] coordinates to NURBS knot span

    u = ...
        (1-xi)*el.u(1)/2 ...
        +(1+xi)*el.u(2)/2;

    v = ...
        (1-eta)*el.v(1)/2 ...
        +(1+eta)*el.v(2)/2;


    %% NURBS basis derivatives in u

    dNu = dersbasisfuns( ...
        el.iu, ...
        u, ...
        el.order(1)-1, ...
        2, ...
        nurbs.knots.U{el.patch});


    %% NURBS basis derivatives in v

    dNv = dersbasisfuns( ...
        el.iv, ...
        v, ...
        el.order(2)-1, ...
        2, ...
        nurbs.knots.V{el.patch});


    %% Rational NURBS surface derivatives

    [~,s] = derRat2DBasisFuns( ...
        dNu, ...
        dNv, ...
        el.order(1), ...
        el.order(2), ...
        el.CP, ...
        2, ...
        2);

end


%% ========================================================================
%  LOCAL FUNCTION:
%  DELTA EVALUATION
%  ========================================================================

function d = ip1_delta(edgeLength,referenceLength,N)

    % Lobatto coordinates along physical edge

    x = ...
        (edgeLength/2)*lobat(N);


    % Normalized maximum Lobatto spacing

    d = ...
        max(diff(x)) ...
        /sqrt(referenceLength*edgeLength);

end


%% ========================================================================
%  LOCAL FUNCTION:
%  SEM OPERATORS
%  ========================================================================

function op = ip1_operators(nx,ny)

    %% Xi direction

    space = struct( ...
        'a',-1, ...
        'b',1, ...
        'N',nx);

    [Fx,Bx] = cheb(space);

    Dx = derivative(space);

    Vx = InnerProduct(space);


    %% Eta direction

    space.N = ny;

    [Fy,By] = cheb(space);

    Dy = derivative(space);

    Vy = InnerProduct(space);


    %% Tensor-product SEM operators

    op.FT = ...
        kron(Fx,Fy);

    op.BT = ...
        kron(Bx,By);

    op.VD = ...
        kron(Vx,Vy);

    op.Q1xi = ...
        kron( ...
        Bx*Dx*Fx, ...
        eye(ny));

    op.Q1eta = ...
        kron( ...
        eye(nx), ...
        By*Dy*Fy);

    op.Q2xi = ...
        kron( ...
        Bx*Dx^2*Fx, ...
        eye(ny));

    op.Q2eta = ...
        kron( ...
        eye(nx), ...
        By*Dy^2*Fy);

end


%% ========================================================================
%  LOCAL FUNCTION:
%  CLOSE FILES OPENED BY LEGACY GEOMETRY READER
%  ========================================================================

function ip1_close_new_files(before)

    handles = ...
        setdiff(openedFiles,before);

    for h = handles(:).'

        fclose(h);

    end

end


%% ========================================================================
%  LOCAL FUNCTION:
%  PLOT NURBS GEOMETRY AND SEM POINTS
%  ========================================================================

function ip1_plot_mesh(nurbs,mesh,nurbsPlotResolution)

    %% --------------------------------------------------------------------
    %  Figure setup
    %  --------------------------------------------------------------------

    tag = 'IP1_SEM_MESH';

    delete( ...
        findall( ...
        groot, ...
        'Type','figure', ...
        'Tag',tag));

    f = figure( ...
        'Name','IP1 - Delta-controlled SEM mesh', ...
        'Tag',tag, ...
        'Color','w', ...
        'Position',[100 100 1350 620]);

    tiledlayout( ...
        f, ...
        1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    %% First axes

    ax1 = nexttile;

    hold(ax1,'on');


    %% Second axes

    ax2 = nexttile;

    hold(ax2,'on');


    %% ====================================================================
    %  FIGURE 1:
    %  NURBS KNOT-SPAN MESH
    %  ====================================================================

    for e = 1:mesh.nel

        xyz = ...
            mesh.boundary{e};


        %% Exact NURBS element boundary

        plot3( ...
            ax1, ...
            xyz(:,1), ...
            xyz(:,2), ...
            xyz(:,3), ...
            '-', ...
            'Color','k', ...
            'LineWidth',1.2);


        %% Element number

        center = ...
            mean( ...
            mesh.nodes(mesh.nodeConn{e},:), ...
            1);

        text( ...
            ax1, ...
            center(1), ...
            center(2), ...
            center(3), ...
            sprintf('%d',e), ...
            'FontSize',12, ...
            'HorizontalAlignment','center');

    end


    %% Figure 1 title

    title( ...
        ax1, ...
        sprintf( ...
        'NURBS coarse mesh: %d elements', ...
        mesh.nel));


    %% ====================================================================
    %  FIGURE 2:
    %  EXACT NURBS SURFACE + SEM POINTS
    %  ====================================================================

    %% Dense plotting coordinates

    xiPlot = ...
        linspace( ...
        -1,1,nurbsPlotResolution);

    etaPlot = ...
        linspace( ...
        -1,1,nurbsPlotResolution);


    %% --------------------------------------------------------------------
    %  Loop over NURBS knot-span elements
    %  --------------------------------------------------------------------

    for e = 1:mesh.nel

        el = ...
            mesh.geometryElements(e);


        %% ---------------------------------------------------------------
        %  Dense exact NURBS surface
        %  ---------------------------------------------------------------

        Xn = zeros( ...
            nurbsPlotResolution, ...
            nurbsPlotResolution);

        Yn = zeros( ...
            nurbsPlotResolution, ...
            nurbsPlotResolution);

        Zn = zeros( ...
            nurbsPlotResolution, ...
            nurbsPlotResolution);


        %% Evaluate exact NURBS surface

        for i = 1:nurbsPlotResolution

            for j = 1:nurbsPlotResolution

                s = ip1_surface( ...
                    nurbs, ...
                    el, ...
                    xiPlot(i), ...
                    etaPlot(j));

                Xn(j,i) = ...
                    s(1,1,1);

                Yn(j,i) = ...
                    s(2,1,1);

                Zn(j,i) = ...
                    s(3,1,1);

            end

        end


        %% ---------------------------------------------------------------
        %  Plot exact NURBS surface
        %  ---------------------------------------------------------------

        surf( ...
            ax2, ...
            Xn, ...
            Yn, ...
            Zn, ...
            'FaceColor',[0.8 1.0 0.8], ...
            'EdgeColor','none', ...
            'FaceAlpha',1.0);


        %% ---------------------------------------------------------------
        %  Exact curved NURBS element boundaries
        %  ---------------------------------------------------------------

        t = ...
            linspace( ...
            -1,1,nurbsPlotResolution).';

        edgeParams = [ ...
            t,               -ones(size(t));
            ones(size(t)),    t;
            flipud(t),        ones(size(t));
            -ones(size(t)),   flipud(t)];

        boundaryXYZ = ...
            zeros(size(edgeParams,1),3);


        %% Evaluate boundary exactly from NURBS

        for a = 1:size(edgeParams,1)

            s = ip1_surface( ...
                nurbs, ...
                el, ...
                edgeParams(a,1), ...
                edgeParams(a,2));

            boundaryXYZ(a,:) = ...
                s(:,1,1).';

        end


        %% Plot curved NURBS boundaries

        plot3( ...
            ax2, ...
            boundaryXYZ(:,1), ...
            boundaryXYZ(:,2), ...
            boundaryXYZ(:,3), ...
            '-', ...
            'Color','k', ...
            'LineWidth',0.8);


        %% ---------------------------------------------------------------
        %  SEM / Lobatto points
        %  ---------------------------------------------------------------

        xyzSEM = ...
            mesh.nodes( ...
            mesh.nodeConn{e},:);


        %% Plot SEM points on top of exact NURBS geometry

        plot3( ...
            ax2, ...
            xyzSEM(:,1), ...
            xyzSEM(:,2), ...
            xyzSEM(:,3), ...
            '.', ...
            'Color','b', ...
            'MarkerSize',15);

    end


    %% ====================================================================
    %  FIGURE 2 TITLE
    %  ====================================================================

    if mesh.samePolynomialDegree

        % In same-degree mode, all elements are required to have one
        % common polynomial order.

        elementOrders = unique(mesh.polynomialOrder(:));

        assert(numel(elementOrders)==1, ...
            ['Expected one common element order, but multiple ' ...
             'polynomial orders were detected.']);

        elementOrder = elementOrders(1);

        % title( ...
        %     ax2, ...
        %     sprintf( ...
        %     'SEM-BEM points: \\delta = %.3g, Element Order = %d', ...
        %     mesh.delta0, ...
        %     elementOrder));

    else

        % Free mode can contain different orders in xi and eta and/or
        % between elements.

        minOrder = min(mesh.polynomialOrder(:));
        maxOrder = max(mesh.polynomialOrder(:));

        if minOrder == maxOrder

            title( ...
                ax2, ...
                sprintf( ...
                'SEM-BEM points: \\delta = %.3g, Element Order = %d', ...
                mesh.delta0, ...
                minOrder));

        else

            title( ...
                ax2, ...
                sprintf( ...
                ['SEM-BEM points: \\delta = %.3g, ' ...
                 'Element Orders = %d-%d'], ...
                mesh.delta0, ...
                minOrder, ...
                maxOrder));

        end

    end


    %% ====================================================================
    %  COMMON AXES SETTINGS
    %  ====================================================================

    for ax = [ax1,ax2]

        axis(ax,'equal');

        axis(ax,'off');

        if max(mesh.nodes(:,3)) ...
                - min(mesh.nodes(:,3)) < 1e-12

            view(ax,2);

        else

            view(ax,3);

        end

    end

end