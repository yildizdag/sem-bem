clear
close all
clc

tic

%% ------------------------------------------------------------------------
% ------------------- Input Parameters ------------------------------------
% -------------------------------------------------------------------------

load SquarePanel_wCenterHole
elements = conn; 
clear conn

% Plotting the geometry using the nodes coordinates
figure(100)
locA = (1:1:max(max(elements)))';
for i=1:size(elements,1)
    plot(nodes(elements(i,:),1), nodes(elements(i,:),2), 'o')
    hold on
end
text(nodes(:,1),nodes(:,2),num2str(locA),'FontSize',15)
hold off
axis equal

% Material properties (Aluminum)
rho_plate = 2700;
E_plate = 69e9;
pois_plate = 0.33;

% Dimensions, h_plate is thickness of the plate. 
h_plate = 0.005;                  % thickness of the host structure (m)


% Boundary conditions
% Rows 1-2 : define the edge points (i.e. the corner nodes in Rhino data)
% Rows 3-7 : define the DOFs
%   - Row 3: to constrain u DOF
%   - Row 4: to constrain v DOF
%   - Row 5: to constrain w DOF
%   - Row 6: to constrain phi_x DOF
%   - Row 7: to constrain phi_y DOF

% Unconstrained Panel:
% BC_edgeinfo = [ ];  

% Fully-clamped Panel:
% BC_edgeinfo = [ 1 76 1 1 1 1 1; 
%                 1 5 1 1 1 1 1; 
%                 5 80 1 1 1 1 1; 
%                 76 80 1 1 1 1 1];

% Simply-supported Panel:
BC_edgeinfo = [ 1 76 1 1 1 0 1; 
                1 5 1 1 1 1 0; 
                5 80 1 1 1 0 1; 
                76 80 1 1 1 1 0];

%% ------------------------------------------------------------------------
% ------------------- Element Sampling ------------------------------------
% -------------------------------------------------------------------------
                            
[indAR, elementpoints, polynums] = element_sampling(elements, nodes, BC_edgeinfo);
% 1) indAR = [indA indR indB indC indD]
%       - indA: DOF information (elements shows the corresponding dofs (i.e., 1-6) 
%       considering the number of sampling of points
%       - indR: shared/unshared sampling points
%       value is 1 if the sampling point is shared by different elements
%       - indB: DOFs (1:N) repeated 6 times
%       - indC: indices for DOFS to be fixed
%       - indD: tangent information for the DOFs to be simply supported 
% 2) elementpoints: index positions of the element (sampling) points
% 3) polynums = [ polynum_xi, polynum_eta]: 
%       polynom numbers for each element along x and y directions

%% ------------------------------------------------------------------------
% -------------- Assembly -------------------------------------------------
% -------------------------------------------------------------------------

% Initializing the system matrices as sparse matrices 
Ka = sparse(size(indAR,1),size(indAR,1));
Ma = sparse(size(indAR,1),size(indAR,1));

% Initializing positions (x, y, z coordinates) for all sampling points in
% the assembly
posn = zeros(size(indAR,1),3);

for di1 = 1:size(elements,1)
    
    % locs: location array (25,2)
    % xlocalnow, ylocalnow: unit vector in local x- and y-directions
    % indelm: indices of sampling points in the assembly matrices
    % Tnow2: transformation matrix
    [locs, xlocalnow, ylocalnow, indelm, Tnow2] = ...
        element_prepare(elements(di1,:),nodes, elementpoints(di1,:), indAR);
    
    % Kelm, Melm: system matrices in local domain
    % xelm, yelm: x and y location in matrix form
    [xelm,yelm,Kelm,Melm] = ...
        Mass_and_Stiffness_Element2(rho_plate,E_plate,pois_plate,2,2,...
        h_plate, polynums(di1,:),locs);
    
    % Distance between the local-global system (needs to be determined ????)
    xyz0now = [0 0 0];
    posnelm = [xelm(:) yelm(:) zeros(length(xelm(:)),1)];
    zlocalnow = cross(xlocalnow,ylocalnow);    
    Tnow = [xlocalnow; ylocalnow; zlocalnow];
    
    % Applying transformation to the coordinates in the local coordinates
    posnnow = ( Tnow \ posnelm')';
    % Adding the distance between the local and global coordinate systems
    posnnow(:,1) = posnnow(:,1)+xyz0now(1);
    posnnow(:,2) = posnnow(:,2)+xyz0now(2);
    posnnow(:,3) = posnnow(:,3)+xyz0now(3);
    % Rounding the coordinates to 12 digits right of the decimal point 
    % (since there may be small numerical differences in the calculated coordinates)
    posnnow = round(posnnow,12);

    % Extending the position vector to all 6-DOFs
    posn(indelm,:) = repmat(posnnow,6,1);

    % Plotting the sampling point on the element in the global coordinate system
    figure(101)
    plot(posnelm(:,1),posnelm(:,2), 'ro')
    hold on 
    locB = (1:1:size(posnelm,1))';
    % hold off
    text(posnelm(:,1),posnelm(:,2), num2str(locB),'FontSize',10)
    axis equal
    
    % Assembly of the system matrices
    Kelmnow = [Kelm     zeros(size(Kelm,1),size(Kelm,1)/5);
               zeros(size(Kelm,1)/5,size(Kelm,1)*6/5)];
    Ka(indelm,indelm) = Ka(indelm,indelm) + Tnow2 \ Kelmnow * Tnow2;
    
    Melmnow = [Melm zeros(size(Kelm,1),size(Kelm,1)/5);
        zeros(size(Kelm,1)/5,size(Kelm,1)*6/5)];
    Ma(indelm,indelm) = Ma(indelm,indelm) + Tnow2 \ Melmnow * Tnow2;
    
end

% Removal of unused indices. This is mainly due to lack of one dof (phi_z) 
indF = find(sum(abs(Ma))==0);
Ka(indF,:) = [];
Ka(:,indF) = [];
Ma(indF,:) = [];
Ma(:,indF) = [];
indAR(indF,:) = [];
posn(indF,:) = [];

% Applying boundary conditions: 
%    deleting the indices constructed in indC = indAR(:,4)
Ka(indAR(:,4)==1,:) = [];
Ka(:,indAR(:,4)==1) = [];
Ma(indAR(:,4)==1,:) = [];
Ma(:,indAR(:,4)==1) = [];
posn(indAR(:,4)==1,:) = [];
indAR(indAR(:,4)==1,:) = [];

disp(['Assembly: ' num2str(round(toc,1)) ' s'])

%% ------------------------------------------------------------------------
% -------------- Eigenvalue Solution --------------------------------------
% -------------------------------------------------------------------------

tic
sigma = 0.01; 
[eigVec,eigVal,flag] = eigs(Ka,Ma,20,sigma);
[wns,loc] = sort(real(sqrt(diag(eigVal)-sigma)));
eigVec = eigVec(:,loc);
disp((wns(1:20))/2/pi)
disp(['Solution time: ' num2str(round(toc,1)) ' s'])
