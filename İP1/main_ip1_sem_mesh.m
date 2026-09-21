%% IP-1: NURBS geometry -> delta-controlled SEM mesh
% Open this file in MATLAB and press Run. No external working-folder setup.
% Geometry evaluation follows sem_core/sem2Dmesh.m.
% Delta selection follows fromBekir_May25/element_sampling.m:
%   - reference length = global x extent;
%   - edge length = straight distance between its corner points;
%   - opposite/shared edges have a common sample count;
%   - each group uses its longest edge;
%   - increase the sample count by one until delta <= delta0.
% The geometric NURBS degree is not changed by delta.
% This example generates a mesh only; it does not call the SEM/BEM solver.
clear; clc;

%% User settings
delta0 = 0.3;
minPolynomialOrder = 1;  % Degree p: p+1 points/polynomials (Bekir: polynow=5)
geometryFolder = 'sem_core/geometry';
fileName = 'plate_cutout_';
numPatch = 4;
shell_dof = 6;
mergeTolerance = 1e-5;   % Same relative coordinate tolerance as sem2Dmesh
showFigures = true;

%% Run (paths are resolved relative to this script, not MATLAB's pwd)
scriptFolder = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptFolder);
geometryPrefix = fullfile(repoRoot,geometryFolder,fileName);
[Nurbs2D,sem2D,elementSummary,edgeSummary] = ip1_build_mesh( ...
    repoRoot,geometryPrefix,numPatch,delta0,minPolynomialOrder, ...
    shell_dof,mergeTolerance);

fprintf('\nGeometry: %s (%d patches)\n',fileName,numPatch);
fprintf('delta0 = %.6g; minimum degree = %d (%d points)\n', ...
    delta0,minPolynomialOrder,minPolynomialOrder+1);
fprintf('Reference x length = %.6g\n',sem2D.referenceLength);
fprintf('%d elements, %d edge groups, %d unique SEM nodes, %d DOFs\n', ...
    sem2D.nel,max(sem2D.edgesgroup),sem2D.non,sem2D.dof);
fprintf('Maximum selected delta = %.6g; minimum sampled surface J = %.6g\n', ...
    max(sem2D.edgeDelta),min(cellfun(@min,sem2D.J)));
disp(elementSummary);
fprintf(['Workspace: Nurbs2D, sem2D, elementSummary, edgeSummary.\n' ...
    'sem2D.polynums(e,:) = [Nxi Neta]; degree = polynums-1.\n' ...
    'conn{e}, nodeConn{e}, J{e}, R{e}, operators{e} are element-wise cells.\n' ...
    'These variable-size arrays are not a drop-in input to global2D.m.\n']);
if showFigures
    ip1_plot_mesh(sem2D);
end

%% Local functions
function [nurbs,mesh,summary,edgeSummary] = ip1_build_mesh( ...
        root,prefix,np,delta0,minDegree,ndof,tol)
    validateattributes(delta0,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(minDegree,{'numeric'},{'scalar','integer','>=',1});
    validateattributes(np,{'numeric'},{'scalar','integer','positive'});
    validateattributes(ndof,{'numeric'},{'scalar','integer','positive'});
    validateattributes(tol,{'numeric'},{'scalar','real','finite','positive'});
    oldPath = path;
    pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
    addpath(fullfile(root,'sem_core'),'-begin');
    for k=1:np
        assert(isfile([prefix num2str(k)]),'Missing geometry file: %s%d',prefix,k);
    end
    % The legacy reader does not close its files. Close only its new handles.
    handlesBefore = openedFiles;
    fileCleanup = onCleanup(@() ip1_close_new_files(handlesBefore)); %#ok<NASGU>
    nurbs = iga2Dmesh(prefix,np,1);
    ip1_close_new_files(handlesBefore);

    %% Original knot-span elements and their physical corners
    ne = sum(cellfun(@double,nurbs.nel));
    elements = repmat(struct('patch',0,'localElement',0,'iu',0,'iv',0, ...
        'u',[],'v',[],'CP',[],'order',[]),ne,1);
    corners = zeros(4*ne,3);
    boundary = cell(ne,1);
    e=0;
    for k=1:np
        for localElement=1:nurbs.nel{k}
            e=e+1;
            iu=nurbs.INC{k}(nurbs.IEN{k}(1,localElement),1);
            iv=nurbs.INC{k}(nurbs.IEN{k}(1,localElement),2);
            u=nurbs.knots.U{k}(iu:iu+1);
            v=nurbs.knots.V{k}(iv:iv+1);
            assert(u(2)>u(1) && v(2)>v(1), ...
                'Zero/negative knot span in patch %d, element %d.',k,localElement);
            order=nurbs.order{k};
            CP=nurbs.cPoints{k}(:,iu-order(1)+1:iu,iv-order(2)+1:iv);
            elements(e)=struct('patch',k,'localElement',localElement, ...
                'iu',iu,'iv',iv,'u',u,'v',v,'CP',CP,'order',order);
            cornerParams=[-1 -1;1 -1;1 1;-1 1];
            for a=1:4
                s=ip1_surface(nurbs,elements(e),cornerParams(a,1),cornerParams(a,2));
                corners(4*(e-1)+a,:)=s(:,1,1).';
            end
            % 5x5 geometry points reproduce the reference code's input-node
            % x-extent convention, independently of the selected SEM count.
            [x,y]=ndgrid(linspace(-1,1,5));
            if e==1, geometryNodes=zeros(25*ne,3); end
            for a=1:25
                s=ip1_surface(nurbs,elements(e),x(a),y(a));
                geometryNodes(25*(e-1)+a,:)=s(:,1,1).';
            end
            t=linspace(-1,1,31).';
            edgeParams=[t,-ones(size(t));ones(size(t)),t; ...
                flipud(t),ones(size(t));-ones(size(t)),flipud(t)];
            boundary{e}=zeros(size(edgeParams,1),3);
            for a=1:size(edgeParams,1)
                s=ip1_surface(nurbs,elements(e),edgeParams(a,1),edgeParams(a,2));
                boundary{e}(a,:)=s(:,1,1).';
            end
        end
    end
    Lref=max(geometryNodes(:,1))-min(geometryNodes(:,1));
    assert(isfinite(Lref) && Lref>0, ...
        'Bekir delta rule requires a positive global x extent.');
    [cornerNodes,~,cornerIds]=uniquetol(corners,tol,'ByRows',true);
    elementCorners=reshape(cornerIds,4,ne).';
    % Local edges follow increasing xi/eta, as in element_sampling.m.
    edgeOrder=[1 2;2 3;4 3;1 4];
    rawEdges=zeros(4*ne,2);
    for a=1:4
        rawEdges((a-1)*ne+(1:ne),:)=elementCorners(:,edgeOrder(a,:));
    end
    [edges,~,edgeIds]=unique(sort(rawEdges,2),'rows');
    elementEdges=reshape(edgeIds,ne,4);
    reversed=reshape(rawEdges(:,1)>rawEdges(:,2),ne,4);
    % Connected components express the transitive opposite-edge equivalence
    % used by Bekir's edge grouping (also for long multi-patch chains).
    pairs=[elementEdges(:,[1 3]);elementEdges(:,[2 4])];
    groups=conncomp(graph(pairs(:,1),pairs(:,2),[],size(edges,1))).';
    edgeLength=vecnorm(cornerNodes(edges(:,2),:)-cornerNodes(edges(:,1),:),2,2);
    assert(all(isfinite(edgeLength) & edgeLength>0),'Degenerate coarse edge.');
    groupLength=accumarray(groups,edgeLength,[],@max);
    groupN=zeros(size(groupLength)); groupDelta=groupN;
    for g=1:numel(groupLength)
        N=minDegree+1;
        d=ip1_delta(groupLength(g),Lref,N);
        while d>delta0
            N=N+1;
            assert(N<=10000,'delta0 requires more than 10000 points per edge.');
            d=ip1_delta(groupLength(g),Lref,N);
        end
        groupN(g)=N; groupDelta(g)=d;
    end
    edgeN=groupN(groups);
    polynums=[edgeN(elementEdges(:,1)),edgeN(elementEdges(:,2))];

    %% Evaluate NURBS at the selected Lobatto points (sem2Dmesh ordering)
    mesh=struct('nel',ne,'shell_dof',ndof,'delta0',delta0, ...
        'minPolynomialOrder',minDegree,'referenceLength',Lref, ...
        'polynums',polynums,'polynomialOrder',polynums-1, ...
        'cornerNodes',cornerNodes,'elementCorners',elementCorners, ...
        'edges',edges,'elementedges',[elementEdges,reversed], ...
        'edgesgroup',groups,'edgeLength',edgeLength, ...
        'edgeGroupLength',groupLength(groups),'edgePolynums',edgeN, ...
        'edgeDelta',groupDelta(groups));
    mesh.N=polynums;  % Explicitly [Nxi Neta] per element, not a scalar.
    mesh.boundary=boundary;
    mesh.patch=[elements.patch].'; mesh.localElement=[elements.localElement].';
    offsets=[0;cumsum(prod(polynums,2))];
    nodeData=zeros(offsets(end),3);
    fields={'nodeConn','conn','J','Jmat','InvJmat','Kappa','R','operators','xi','eta'};
    for f=1:numel(fields),mesh.(fields{f})=cell(ne,1);end
    for e=1:ne
        nx=polynums(e,1);ny=polynums(e,2);nn=nx*ny;
        xi=lobat(nx);eta=lobat(ny);
        mesh.xi{e}=xi;mesh.eta{e}=eta;
        mesh.J{e}=zeros(nn,1);mesh.Jmat{e}=zeros(3,2,nn);
        mesh.InvJmat{e}=zeros(2,2,nn);mesh.Kappa{e}=zeros(2,nn);
        mesh.R{e}=zeros(3,3,nn);
        du=diff(elements(e).u)/2;dv=diff(elements(e).v)/2;
        a=0;
        for i=1:nx
            for j=1:ny
                a=a+1;s=ip1_surface(nurbs,elements(e),xi(i),eta(j));
                A1=s(:,2,1);A2=s(:,1,2);A=[A1,A2];metric=A.'*A;
                normal=cross(A1,A2);area=norm(normal);
                assert(isfinite(area) && area>0 && rcond(metric)>eps, ...
                    'Degenerate surface: element %d, sample %d.',e,a);
                normal=normal/area;t1=A1/norm(A1);t2=cross(normal,t1);t2=t2/norm(t2);
                Ac=A/metric;
                curvature=[dot(s(:,3,1),normal),dot(s(:,2,2),normal); ...
                    dot(s(:,2,2),normal),dot(s(:,1,3),normal)];
                nodeData(offsets(e)+a,:)=s(:,1,1).';
                mesh.J{e}(a)=area*du*dv;
                mesh.Jmat{e}(:,:,a)=A;
                mesh.InvJmat{e}(:,:,a)=[dot(t1,Ac(:,1))/du,dot(t1,Ac(:,2))/dv; ...
                    dot(t2,Ac(:,1))/du,dot(t2,Ac(:,2))/dv];
                mesh.Kappa{e}(:,a)=abs(eig(metric\curvature));
                mesh.R{e}(:,:,a)=[t1,t2,normal];
            end
        end
        mesh.operators{e}=ip1_operators(nx,ny);
    end
    % Preserve the existing mesh's coordinate-based global node merging.
    [mesh.nodes,~,ids]=uniquetol(nodeData,tol,'ByRows',true);
    for e=1:ne
        nc=ids(offsets(e)+(1:prod(polynums(e,:)))).';
        mesh.nodeConn{e}=nc;
        conn=zeros(1,ndof*numel(nc));
        for d=1:ndof,conn(d:ndof:end)=ndof*nc-(ndof-d);end
        mesh.conn{e}=conn;
    end
    mesh.non=size(mesh.nodes,1);mesh.dof=ndof*mesh.non;
    % Verify that edge-group compatibility also holds geometrically.
    edgeNodes=cell(size(edges,1),1);
    for e=1:ne
        nx=polynums(e,1);ny=polynums(e,2);
        localEdges={1:ny:(nx-1)*ny+1,(nx-1)*ny+(1:ny),ny:ny:nx*ny,1:ny};
        for a=1:4
            current=mesh.nodeConn{e}(localEdges{a});
            if reversed(e,a),current=fliplr(current);end
            id=elementEdges(e,a);
            if isempty(edgeNodes{id})
                edgeNodes{id}=current;
            else
                assert(isequal(edgeNodes{id},current), ...
                    ['Shared edge %d has incompatible sampled coordinates. ' ...
                     'Matching endpoints alone do not guarantee an identical curve/parameterization.'],id);
            end
        end
    end
    mesh.edgeNodes=edgeNodes;
    deltaXi=mesh.edgeDelta(elementEdges(:,1));deltaEta=mesh.edgeDelta(elementEdges(:,2));
    summary=table((1:ne).',mesh.patch,mesh.localElement,polynums(:,1),polynums(:,2), ...
        polynums(:,1)-1,polynums(:,2)-1,deltaXi,deltaEta,cellfun(@min,mesh.J), ...
        'VariableNames',{'Element','Patch','PatchElement','Nxi','Neta', ...
        'DegreeXi','DegreeEta','DeltaXi','DeltaEta','MinSurfaceJ'});
    edgeSummary=table((1:size(edges,1)).',groups,edgeLength,groupLength(groups), ...
        edgeN,mesh.edgeDelta,'VariableNames', ...
        {'Edge','Group','ChordLength','GroupMaxLength','Points','Delta'});
end

function s=ip1_surface(nurbs,el,xi,eta)
    u=(1-xi)*el.u(1)/2+(1+xi)*el.u(2)/2;
    v=(1-eta)*el.v(1)/2+(1+eta)*el.v(2)/2;
    dNu=dersbasisfuns(el.iu,u,el.order(1)-1,2,nurbs.knots.U{el.patch});
    dNv=dersbasisfuns(el.iv,v,el.order(2)-1,2,nurbs.knots.V{el.patch});
    [~,s]=derRat2DBasisFuns(dNu,dNv,el.order(1),el.order(2),el.CP,2,2);
end

function d=ip1_delta(edgeLength,referenceLength,N)
    % Identical physical points to Discretization(edgeLength,N,'xi').
    x=(edgeLength/2)*lobat(N);
    d=max(diff(x))/sqrt(referenceLength*edgeLength);
end

function op=ip1_operators(nx,ny)
    space=struct('a',-1,'b',1,'N',nx);
    [Fx,Bx]=cheb(space);Dx=derivative(space);Vx=InnerProduct(space);
    space.N=ny;
    [Fy,By]=cheb(space);Dy=derivative(space);Vy=InnerProduct(space);
    op.FT=kron(Fx,Fy);op.BT=kron(Bx,By);op.VD=kron(Vx,Vy);
    op.Q1xi=kron(Bx*Dx*Fx,eye(ny));op.Q1eta=kron(eye(nx),By*Dy*Fy);
    op.Q2xi=kron(Bx*Dx^2*Fx,eye(ny));op.Q2eta=kron(eye(nx),By*Dy^2*Fy);
end

function ip1_close_new_files(before)
    handles=setdiff(openedFiles,before);
    for h=handles(:).',fclose(h);end
end

function ip1_plot_mesh(mesh)
    % Only replace this example's figure; keep the user's other figures.
    tag='IP1_SEM_MESH';delete(findall(groot,'Type','figure','Tag',tag));
    f=figure('Name','IP1 - Delta-controlled SEM mesh','Tag',tag,'Color','w', ...
        'Position',[100 100 1350 620]);
    tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax1=nexttile;hold(ax1,'on');ax2=nexttile;hold(ax2,'on');
    colors=lines(max(mesh.patch));
    for e=1:mesh.nel
        xyz=mesh.boundary{e};c=colors(mesh.patch(e),:);
        plot3(ax1,xyz(:,1),xyz(:,2),xyz(:,3),'-','Color',c,'LineWidth',1.2);
        center=mean(mesh.nodes(mesh.nodeConn{e},:),1);
        text(ax1,center(1),center(2),center(3),sprintf('%d',e),'FontSize',9);
        nx=mesh.polynums(e,1);ny=mesh.polynums(e,2);
        xyz=mesh.nodes(mesh.nodeConn{e},:);
        X=reshape(xyz(:,1),ny,nx);Y=reshape(xyz(:,2),ny,nx);Z=reshape(xyz(:,3),ny,nx);
        surf(ax2,X,Y,Z,'FaceColor','none','EdgeColor',c,'LineWidth',0.5);
        plot3(ax2,xyz(:,1),xyz(:,2),xyz(:,3),'.','Color',c,'MarkerSize',8);
    end
    title(ax1,sprintf('NURBS coarse mesh: %d elements',mesh.nel));
    title(ax2,sprintf('SEM: delta = %.3g, Poly. Degree = %d, %d nodes', ...
        mesh.delta0,mesh.minPolynomialOrder,mesh.non));
    for ax=[ax1,ax2]
        axis(ax,'equal');grid(ax,'on');xlabel(ax,'x');ylabel(ax,'y');zlabel(ax,'z');
        if max(mesh.nodes(:,3))-min(mesh.nodes(:,3))<1e-12,view(ax,2);else,view(ax,3);end
    end
end
