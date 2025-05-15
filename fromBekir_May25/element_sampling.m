function [sem_mesh] = element_sampling(sem_mesh)
% ELEMENT_SAMPLING Generates high-order sampling point data for SEM
%
% This function constructs the data structures required for spectral element
% method (SEM) analysis based on a given mesh consisting of high-order elements.
% Each element is assumed to contain a fixed number of sampling points defined
% according to the discretization scheme shown in Fig. 2 of 
% https://doi.org/10.1016/j.tws.2024.111636. 
%
% INPUT:
%   sem_mesh - Struct containing:
%       .elements : Element connectivity matrix (each row corresponds to one 
%                   25-node high-order quadrilateral element)
%       .nodes    : Nodal coordinates (Nx3) defining the geometry of the mesh
%
% OUTPUT:
%   sem_mesh - Updated struct including:
%       .elements         - (unchanged) Input element connectivity
%       .nodes            - (unchanged) Input nodal coordinates
%       .elementedgenodes - Corner node indices for each element (Rhino-defined)
%       .edges            - Global unique edge list (with duplicates removed)
%       .elementedges     - Edge indices and orientation flags for each element
%       .edgesgroup       - Grouping of equivalent edges (for polynomial assignment)
%       .polynums         - Polynomial orders in xi and eta directions per element
%       .nodepoints       - Global IDs of corner sampling points (only assigned once)
%       .edgepoints       - Global IDs of edge sampling points (excluding corners)
%       .elementpoints    - Global sampling point indices for each element
%       .group_ind        - Group identifier per sampling point (corner/edge/interior)
%       .ind_ALL          - Global DOF list: [DOF_type, shared_flag, point_index]
%
% The function:
%   - Detects unique edges from corner nodes
%   - Groups edges that are opposite in elements (for polynomial consistency)
%   - Assigns polynomial orders per edge group using a delta-based convergence criterion
%   - Generates and numbers sampling points uniquely across the mesh
%   - Classifies sampling points as node, edge, or interior
%   - Identifies shared DOFs across elements for global system assembly
%--------------------------------------------------------------------------



% calculating the length of the plate along x-direction (may need to change in the future)
Lx_plate = max(sem_mesh.nodes(:,1)) - min(sem_mesh.nodes(:,1));

%% Corner points of each element ------------------------------------------
% The 25 points for each element are defined increasing along x-axis first
% then along y-axis (check Fig.2 in https://doi.org/10.1016/j.tws.2024.111636)
%
%         21 o---------o 25
%            |        /
%            |       /
%          1 o------o 5

% Four corner nodes of each element
sem_mesh.elementedgenodes = sem_mesh.elements(:,[1 5 25 21]);

%% Defining all edges in the assembly ------------------------------------- 

% Edges of elements are created by connecting corner points
edges_raw = [sem_mesh.elementedgenodes(:,[1 2]); sem_mesh.elementedgenodes(:,[2 3]); ...
             sem_mesh.elementedgenodes(:,[3 4]); sem_mesh.elementedgenodes(:,[4 1])];
% There may be repeated edges, so the edges are sorted and repeated ones
% are eliminated
sem_mesh.edges = unique(sort(edges_raw,2),'rows');

%% Defining element edges -------------------------------------------------

% Defining the edge directions 
%        4 o----->---------o 3
%          ↑              ↗  
%          |            ↗    
%        1 o----->-----o 2  
edgeorder = [1 2; 2 3; 4 3; 1 4];
% Initializing elementedge matrix:
% First 4 columns include edge numbers of elements
% Last 4 columns include direction information 
sem_mesh.elementedges = zeros(size(sem_mesh.elementedgenodes,1),8);
for di1 = 1:size(sem_mesh.elementedgenodes,1)
    for di2 = 1:4
        aa = sem_mesh.elementedgenodes(di1,edgeorder(di2,:));
        if aa(1) > aa(2)
            aa = aa([2 1]);
            sem_mesh.elementedges(di1,di2+4) = 1;
        end
        sem_mesh.elementedges(di1,di2) = ...
            find((sem_mesh.edges(:,1)==aa(1)).*(sem_mesh.edges(:,2)==aa(2)));
    end 
end

%% Grouping edges with equivalent polynomail numbers (sampling points) ----

% Setting the opposite edges in each element for the whole assembly
%        4 o-----e7-----o 3
%          |            |  
%         e2            e10
%          |            |    
%        1 o-----e3-----o 2  
%Example: e3 and e7 are opposite and e2 and e10 are opposite to each other
edgeseqv = zeros(size(sem_mesh.edges,1),size(sem_mesh.edges,1));
for di1 = 1:size(sem_mesh.elementedgenodes,1)
    edgeseqv(sem_mesh.elementedges(di1,[1 3]),sem_mesh.elementedges(di1,[1 3])) = 1;
    edgeseqv(sem_mesh.elementedges(di1,[2 4]),sem_mesh.elementedges(di1,[2 4])) = 1;    
end

for di1 = 1:size(edgeseqv,1)-1
    for di2 = di1+1:size(edgeseqv,1)
        if sum(sum(edgeseqv([di1 di2],:),1)==2) > 0
            edgeseqv(di1,di2) = 1;
            edgeseqv(di2,di1) = 1;
        end
    end
end

% Determining polynum equivalent edge groups
sem_mesh.edgesgroup = zeros(size(sem_mesh.edges,1),1);
for di1 = 1:size(sem_mesh.edges,1)                                                               
    aa = find(edgeseqv(di1,1:di1-1)==1);
    if ~isempty(aa)
        sem_mesh.edgesgroup(di1) = sem_mesh.edgesgroup(aa(1));
    else
        sem_mesh.edgesgroup(di1) = max(sem_mesh.edgesgroup)+1;
    end
end


%% Calculating the edge lengths in each group -----------------------------

% To calculate the edge length we used the coordinates of the edge corner,
% and calculate the length of a straight line between these corner points.
% However it should be modified to calculate the length of the curve such that 
% the actual lengths of the edges will be found.

% Lengths of each edge
edgeslength = sqrt((sem_mesh.nodes(sem_mesh.edges(:,1),1) - sem_mesh.nodes(sem_mesh.edges(:,2),1)).^2+...
    (sem_mesh.nodes(sem_mesh.edges(:,1),2) - sem_mesh.nodes(sem_mesh.edges(:,2),2)).^2+...
    (sem_mesh.nodes(sem_mesh.edges(:,1),3) - sem_mesh.nodes(sem_mesh.edges(:,2),3)).^2);

% Assigning edge lengths based on the edge groups
edgeslength1 = zeros(size(edgeslength));
for di1 = 1:max(sem_mesh.edgesgroup)
    aa = find(sem_mesh.edgesgroup==di1);
    edgeslength1(aa) = max(edgeslength(aa));
end

%% Determining polynums for given a delta value (convergence amount) ------

polynum_gen = zeros(size(sem_mesh.edges,1),1);
for di1 = 1:length(polynum_gen)
    delta0 = 0.25;
    polynow = 5;
    [~,~,~,xelm,~,~,~,~] = Discretization(edgeslength1(di1),polynow,'xi');
    deltanowx = max(xelm(2:end)-xelm(1:end-1))/sqrt(Lx_plate*edgeslength1(di1));
    while deltanowx > delta0
        polynow = polynow+1;
        [~,~,~,xelm,~,~,~,~] = Discretization(edgeslength1(di1),polynow,'xi');
        deltanowx = max(xelm(2:end)-xelm(1:end-1))/sqrt(Lx_plate*edgeslength1(di1));
    end
    polynum_gen(di1) = polynow;
end

% polynums for each element
polynum_xi = polynum_gen(sem_mesh.elementedges(:,1));
polynum_eta = polynum_gen(sem_mesh.elementedges(:,2));

sem_mesh.polynums = [polynum_xi polynum_eta];


%% Grouping sampling points in each element ------------------------------- 
%
% For each node, nodepoints array is checked if a number is assigned before
% if yes, it is used for each edge, 
% edgepoints array is checked if numbers are assigned before
% if yes, they are used
%
% all unnumbered nodes are given numbers. 
%
% all nodes and edges are then assigned the numbers into nodepoints and 
% edgepoints arrays. 
%
% elementpoints array is created to indicate the indices of the 
% sampling points on each element. 


sem_mesh.nodepoints = zeros(size(sem_mesh.nodes,1),1);
sem_mesh.edgepoints = zeros(size(sem_mesh.edges,1),max(polynum_gen));
sem_mesh.elementpoints = zeros(size(sem_mesh.elementedgenodes,1),max(polynum_xi.*polynum_eta));
sem_mesh.group_ind = zeros(size(sem_mesh.elementedgenodes,1),max(polynum_xi.*polynum_eta));

count = 1;
for di1 = 1:size(sem_mesh.elementedgenodes,1)
    
    % For each element, the point indices are grouped into 9 groups: 
    %   - 4 nodes (corner points)
    %   - 4 edges (sampling points on the edges excluding the corner points)
    %   - 1 interior group (all the sampling points inside the edges) 
    %
    %  4   7   3
    %  .-------.
    %  |       |
    %  |       |
    % 8|   9   |6
    %  |       |
    %  |       |
    %  .-------.
    %  1   5   2
    %
    groupindex = zeros(1,polynum_xi(di1)*polynum_eta(di1));
    groupindex(1) = 1;
    groupindex((polynum_xi(di1)-1)*polynum_eta(di1)+1) = 2;
    groupindex(polynum_xi(di1)*polynum_eta(di1)) = 3;
    groupindex(polynum_eta(di1)) = 4;
    groupindex((polynum_eta(di1)+1):(polynum_eta(di1)):((polynum_xi(di1)-2)*polynum_eta(di1)+1)) = 5;
    groupindex(((polynum_xi(di1)-1)*polynum_eta(di1)+2):(polynum_xi(di1)*polynum_eta(di1)-1)) = 6;
    groupindex((polynum_eta(di1)*2):(polynum_eta(di1)):((polynum_xi(di1)-1)*polynum_eta(di1))) = 7;
    groupindex(2:(polynum_eta(di1)-1)) = 8;    
    groupindex(groupindex==0) = 9;

    sem_mesh.group_ind(di1,1:numel(groupindex)) = groupindex;
    
    pointsnow = zeros(1,polynum_xi(di1)*polynum_eta(di1));
    polynumsnow = [polynum_xi(di1) polynum_eta(di1) polynum_xi(di1) polynum_eta(di1)];
    
    % checking nodes and edges if previously numbered
    % For each node (corner points of elements), nodepoints array is
    % checked if a number is assigned before: if yes, it is used 
    % Similarly for each edge edgepoints array is checked if numbers are
    % assigned before: if yes, they are used
    for di2 = 1:4
        if sem_mesh.nodepoints(sem_mesh.elementedgenodes(di1,di2))>0
            pointsnow(groupindex==di2) = sem_mesh.nodepoints(sem_mesh.elementedgenodes(di1,di2));
        end
    end
    for di2 = 1:4
        if sem_mesh.edgepoints(sem_mesh.elementedges(di1,di2),1)>0
            if sem_mesh.elementedges(di1,di2+4) == 0
                pointsnow(groupindex==di2+4) = sem_mesh.edgepoints(sem_mesh.elementedges(di1,di2),1:polynumsnow(di2)-2);
            else
                pointsnow(groupindex==di2+4) = sem_mesh.edgepoints(sem_mesh.elementedges(di1,di2),polynumsnow(di2)-2:-1:1);                
            end
        end
    end
    
    % For the considered elements, unnumbered (sampling) points are numbered
    count1 = count+sum(pointsnow==0);
    pointsnow(pointsnow==0) = count:(count1-1);
    % The corresponding sampling point numbers are assigned to the elements
    sem_mesh.elementpoints(di1,1:length(pointsnow)) = pointsnow;
    count = count1;
    
    % For the whole assembly, (corner points) nodes and edges are numbered
    for di2 = 1:4
        sem_mesh.nodepoints(sem_mesh.elementedgenodes(di1,di2)) = pointsnow(groupindex==di2);
    end
    for di2 = 1:4
        if sem_mesh.elementedges(di1,di2+4) == 0
            sem_mesh.edgepoints(sem_mesh.elementedges(di1,di2),1:polynumsnow(di2)-2) = pointsnow(groupindex==di2+4);
        else
            aa = pointsnow(groupindex==di2+4);
            sem_mesh.edgepoints(sem_mesh.elementedges(di1,di2),1:polynumsnow(di2)-2) = aa(end:-1:1);
        end
    end
    
end

%% DOF information: -------------------------------------------------------
% indA: DOF information (elements shows the corresponding dofs (i.e., 1-6) 
%       considering the number of sampling of points
% indR: shared/unshared sampling points
%       value is 1 if the sampling point is shared by different elements
% indB: DOFs (1:N) repeated 6 times

% finding how many times each sampling point is used
% if a sampling point is shared between only two elements, the value will be 2
% if a sampling point is shared between four elements, the value will be 4
indr = zeros(max(sem_mesh.elementpoints(:)),1);
for di1 = 1:size(sem_mesh.elementedgenodes,1)
    indnow = sem_mesh.elementpoints(di1,1:(polynum_xi(di1)*polynum_eta(di1)));
    indr(indnow) = indr(indnow)+1;    
end
% Assigning the value 1 for shared sampling points
indR = indr>1;
% Extending the shared sampling points to the DOFs 
indR = repmat(indR,6,1);

% Construction the DOF information for sampling points 
% (u->1; v->2; w->3; phi_x->4; phi_y->5; phi_z->6)
indA = [ones(length(indr),1); 2*ones(length(indr),1); ...
      3*ones(length(indr),1); 4*ones(length(indr),1); ...
      5*ones(length(indr),1); 6*ones(length(indr),1) ];

% Sampling point numbering repeated for DOFs (u,v,w,phi_x,phi_y,phi_z)
% [ (1:Nx*Ny)'; (1:Nx*Ny)'; ...; (1:Nx*Ny)']
indB = repmat((1:length(indr))',6,1);

sem_mesh.ind_ALL = [indA indR indB];

end