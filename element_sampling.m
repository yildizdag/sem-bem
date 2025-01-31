function [indA,indB,indR,elementpoints,polynum_xi,polynum_eta] = ...
    element_sampling(elements,nodes)

%% ------------------------------------------------------------------------
% -------------- Element Sampling -----------------------------------------
% -------------------------------------------------------------------------

Lx_plate = max(nodes(:,1))-min(nodes(:,1));

elementedgenodes = elements(:,[1 5 25 21]);

% edges: all edges in the assembly
edges = [elementedgenodes(:,[1 2]); elementedgenodes(:,[2 3]); elementedgenodes(:,[3 4]); 
    elementedgenodes(:,[4 1])];
edges = unique(sort(edges,2),'rows');

% elementedges: edge numbers on elements(1:4), direction (5:8)
edgeorder = [1 2; 2 3; 4 3; 1 4];
elementedges = zeros(size(elementedgenodes,1),8);
for di1 = 1:size(elementedgenodes,1)
    for di2 = 1:4
        aa = elementedgenodes(di1,edgeorder(di2,:));
        if aa(1) > aa(2)
            aa = aa([2 1]);
            elementedges(di1,di2+4) = 1;
        end
        elementedges(di1,di2) = find((edges(:,1)==aa(1)).*(edges(:,2)==aa(2)));
    end 
end

% edges checklist for polynum equivalence
edgeseqv = zeros(size(edges,1),size(edges,1));
for di1 = 1:size(elementedgenodes,1)
    edgeseqv(elementedges(di1,[1 3]),elementedges(di1,[1 3])) = 1;
    edgeseqv(elementedges(di1,[2 4]),elementedges(di1,[2 4])) = 1;    
end

% polynum equivalent edge groups
edgesgroup = zeros(size(edges,1),1);
for di1 = 1:size(edges,1)
    aa = find(edgeseqv(di1,1:di1-1)==1);
    if ~isempty(aa)
        edgesgroup(di1) = edgesgroup(aa(1));
    else
        edgesgroup(di1) = max(edgesgroup)+1;
    end
end

% longest edgelength in each group

edgeslength = sqrt((nodes(edges(:,1),1) - nodes(edges(:,2),1)).^2+...
    (nodes(edges(:,1),2) - nodes(edges(:,2),2)).^2+...
    (nodes(edges(:,1),3) - nodes(edges(:,2),3)).^2);

edgeslength1 = zeros(size(edgeslength));
for di1 = 1:max(edgesgroup)
    aa = find(edgesgroup==di1);
    edgeslength1(aa) = max(edgeslength(aa));
end

% polynum for given delta
polynum_gen = zeros(size(edges,1),1);
for di1 = 1:length(polynum_gen)
delta0 = 0.5;
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
polynum_xi = polynum_gen(elementedges(:,1));
polynum_eta = polynum_gen(elementedges(:,2));



% the arrays nodepoints, edgepoints and elementpoints are used for sampling 
% point equivalence. 

% for each element, the point indices are grouped into 9 groups: 4 nodes,
% 4 edges and 1 interior group. 

%  4   7   3
%  .-------.
%  |       |
%  |       |
% 8|   9   |6
%  |       |
%  |       |
%  .-------.
%  1   5   2

% for each node, nodepoints array is checked if a number is assigned before
% if yes, it is used
% for each edge, edgepoints array is checked if numbers are assigned before
% if yes, they are used

% all unnumbered nodes are given numbers. 

% all nodes and edges are then assigned the numbers into nodepoints and 
% edgepoints arrays. 

% elementpoints array is created to indicate the indices of the 
% sampling points on each element. 


nodepoints = zeros(size(nodes,1),1);
edgepoints = zeros(size(edges,1),max(polynum_gen));
elementpoints = zeros(size(elementedgenodes,1),max(polynum_xi.*polynum_eta));

count = 1;
for di1 = 1:size(elementedgenodes,1)
    
    % grouping
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
    
    pointsnow = zeros(1,polynum_xi(di1)*polynum_eta(di1));
    polynumsnow = [polynum_xi(di1) polynum_eta(di1) polynum_xi(di1) polynum_eta(di1)];
    
    % checking nodes and edges if previously numbered
    for di2 = 1:4
        if nodepoints(elementedgenodes(di1,di2))>0
            pointsnow(groupindex==di2) = nodepoints(elementedgenodes(di1,di2));
        end
    end
    for di2 = 1:4
        if edgepoints(elementedges(di1,di2),1)>0
            if elementedges(di1,di2+4) == 0
                pointsnow(groupindex==di2+4) = edgepoints(elementedges(di1,di2),1:polynumsnow(di2)-2);
            else
                pointsnow(groupindex==di2+4) = edgepoints(elementedges(di1,di2),polynumsnow(di2)-2:-1:1);                
            end
        end
    end
    
    % numbering unnumbered points
    count1 = count+sum(pointsnow==0);
    pointsnow(pointsnow==0) = count:(count1-1);
    elementpoints(di1,1:length(pointsnow)) = pointsnow;
    count = count1;
    
    % numbering nodes and edges
    for di2 = 1:4
        nodepoints(elementedgenodes(di1,di2)) = pointsnow(groupindex==di2);
    end
    for di2 = 1:4
        if elementedges(di1,di2+4) == 0
            edgepoints(elementedges(di1,di2),1:polynumsnow(di2)-2) = pointsnow(groupindex==di2+4);
        else
            aa = pointsnow(groupindex==di2+4);
            edgepoints(elementedges(di1,di2),1:polynumsnow(di2)-2) = aa(end:-1:1);
        end
    end
    
end



indr = zeros(max(elementpoints(:)),1);
for di1 = 1:size(elementedgenodes,1)
    indnow = elementpoints(di1,1:(polynum_xi(di1)*polynum_eta(di1)));
    indr(indnow) = indr(indnow)+1;     
end
indR = indr>1;

indA = [ones(size(indR,1),1); 2*ones(size(indR,1),1); ...
    3*ones(size(indR,1),1); 4*ones(size(indR,1),1); ...
    5*ones(size(indR,1),1); 6*ones(size(indR,1),1) ];

indB = repmat((1:length(indR))',6,1);


% indR: 1 if the sampling point is shared by different elements
indR = repmat(indR,6,1);

