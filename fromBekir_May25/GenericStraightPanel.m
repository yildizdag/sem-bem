clear
close all
clc

tic
profile on

%% ------------------------------------------------------------------------
% ------------------- Geometry Input --------------------------------------
% -------------------------------------------------------------------------

% Load geometry data (SquarePanel_wCenterHole, rhinodata_RectPanelwCutout,
% rhinodata_circularpanel)
load rhinodata_dummy
sem_mesh.elements = conn; 
sem_mesh.nodes = nodes;
clear conn nodes

% Plotting the geometry using the nodes' coordinates
figure(100)
locA = (1:1:max(max(sem_mesh.elements)))';
for i=1:size(sem_mesh.elements,1)
    plot(sem_mesh.nodes(sem_mesh.elements(i,:),1), sem_mesh.nodes(sem_mesh.elements(i,:),2), 'o')
    hold on
end
text(sem_mesh.nodes(:,1), sem_mesh.nodes(:,2), num2str(locA),'FontSize',15)
hold off
axis equal

% Thickness of the host structure (m) 
h_plate = 0.1;

%% ------------------------------------------------------------------------
% ------------------- Material Input --------------------------------------
% -------------------------------------------------------------------------

% Material properties (Aluminum)
rho_plate = 8000;       % density [kg/m3]
E_plate = 210e9;         % Elastic modulus [Pa]
pois_plate = 0.3;      % Poisson's ratio

%% ------------------------------------------------------------------------
% ------------------- Boundary Input --------------------------------------
% -------------------------------------------------------------------------

% Define boundary conditions for edge nodes
% -----------------------------------------
% The BC struct includes:
%   - edgenodes: corner nodes (from Rhino) that define boundary edges
%   - type: boundary condition type for each edge
%            'C' - Clamped
%            'S' - Simply supported
%            'P' - Pinned
%            'F' - Free
%   - fixeddofs: defines which DOFs are fixed:
%       Column 1: u (in-plane x translation)
%       Column 2: v (in-plane y translation)
%       Column 3: w (out-of-plane translation)
%       Column 4: φₓ or φₙ (rotation about y-axis or normal axis)
%       Column 5: φᵧ or φₜ (rotation about x-axis or tangent axis)


% Define the edge corner points (in Rhino data)
% BC.edgenodes = [1 5; 1 21; 5 25; 21 25];
BC.edgenodes = [1 5; 5 9; 9 13; 13 17; 17 21; 21 25; 25 29; 29 33; 33 37; 37 41]; %;  ...
                % 1 165; 165 329; 329 493; 493 657; 657 821; 821 985; 985 1149; 1149 1313; 1313 1477; 1477 1641;  ... 
                % 41 205; 205 369; 369 533; 533 697; 697 861; 861 1025; 1025 1189; 1189 1353; 1353 1517; 1517 1681; ...
                % 1641 1645; 1645 1649; 1649 1653; 1653 1657; 1657 1661; 1665 1669; 1669 1673; 1673 1677; 1677 1681];

% BC.edgenodes = [1 18; 
%                 1 19;
%                 18 97;
%                 19 113; 
%                 97 191;                 
%                 113 192; 
%                 192 209;               
%                 191 209];

% Define boundary condition type for each edge
BC.type = ['S';'S';'S';'S';...
            'S';'S';'S';'S';...
            'S';'S'] ;% ;'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S';...
            % 'S';'S';'S';'S'];
%BC.type = ['C';'C';'C';'S';'S';'S';'S';'S'];

% Map each BC type to its fixed DOF pattern
bc_map = containers.Map({'C', 'S', 'P', 'F'}, {
    [1 1 1 1 1];  % Clamped
    [1 1 1 1 0];  % Simply supported
    [1 1 1 0 0];  % Pinned
    [0 0 0 0 0]   % Free
});

% Assign fixed DOFs based on BC type
BC.fixeddofs = zeros(length(BC.type), 5);
for di1 = 1:length(BC.type)
    BC.fixeddofs(di1,:) = bc_map(BC.type(di1));
end

% Remove rows with no constraints (i.e., 'F' type)
BC.fixeddofs(~any(BC.fixeddofs, 2), :) = [];

%% ------------------------------------------------------------------------
% ------------------- Element Sampling ------------------------------------
% -------------------------------------------------------------------------
                            
% -------------------------------------------------------------------------
% Generate SEM mesh data structure
% This step processes the input nodal coordinates and element connectivity 
% to construct a spectral element mesh (sem_mesh). It:
%   - Identifies unique edges and groups them based on adjacency
%   - Determines polynomial orders for each edge using a convergence criterion
%   - Assigns global sampling point indices (nodes, edges, interior)
%   - Constructs DOF mapping and shared point information
%
% Required inputs:
%   elements : connectivity matrix for high-order elements (25 nodes/element)
%   nodes    : Nx3 coordinate matrix (from Rhino or CAD preprocessor)
% which are included in sem_mesh struct.
%
% Returns:
%   sem_mesh : struct containing all sampling, edge, and DOF-related data
% -------------------------------------------------------------------------

[sem_mesh] = element_sampling(sem_mesh);

%% ------------------------------------------------------------------------
% -------------- Assembly -------------------------------------------------
% -------------------------------------------------------------------------

% Initializing the system matrices as sparse matrices 
Ka = sparse(size(sem_mesh.ind_ALL,1),size(sem_mesh.ind_ALL,1));
Ma = sparse(size(sem_mesh.ind_ALL,1),size(sem_mesh.ind_ALL,1));

% Initializing positions (x-y-z coord.) for all sampling points in the assembly
 sem_mesh.posn = zeros(size(sem_mesh.ind_ALL,1),3);

% Loop over all elements to compute/assemble local mass and stiffness matrices
for di1 = 1:size(sem_mesh.elements,1)

    % Prepare element-specific data:
    % - locs: local sampling point coordinates (25,2)
    % - indelm: global indices of local sampling points
    % - Tnow2: transformation matrix for DOF rotation
    [locs, indelm, Tnow2] = ...
        element_prepare(sem_mesh.elements(di1,:), sem_mesh.nodes, ...
                        sem_mesh.elementpoints(di1,:), sem_mesh.ind_ALL, 0);

    % Compute local stiffness (Kelm) and mass (Melm) matrices
    % xelm and yelm are matrices containing x and y coord. of local points
    [xelm, yelm, Kelm, Melm] = ...
        Mass_and_Stiffness_Element2(rho_plate, E_plate, pois_plate, 2, 2,...
                                    h_plate, sem_mesh.polynums(di1,:), locs);

    % Distance between the local-global system (needs to be determined ????)
    xyz0now = [0 0 0];

    % Convert local 2D coordinates to 3D by padding with zeros in z-direction
    posnelm = [xelm(:) yelm(:) zeros(length(xelm(:)),1)];

    % Transform local coordinates to global frame (rotation skipped here)
    posnnow = (posnelm')';

    % Adding the distance between the local and global coordinate systems
    posnnow(:,1) = posnnow(:,1)+xyz0now(1);
    posnnow(:,2) = posnnow(:,2)+xyz0now(2);
    posnnow(:,3) = posnnow(:,3)+xyz0now(3);
  
    % Round coordinates to 12 decimal places to avoid numerical inconsistencies
    posnnow = round(posnnow,12);

    % Store 3D positions for each DOF (6 DOFs per node)
    sem_mesh.posn(indelm,:) = repmat(posnnow,6,1);

    % % Plotting the sampling point on the element in the global coordinate system
    % figure(101)
    % plot(posnelm(:,1),posnelm(:,2), 'ro')
    % hold on 
    % locB = (1:1:size(posnelm,1))';
    % % hold off
    % text(posnelm(:,1),posnelm(:,2), num2str(locB),'FontSize',10)
    % axis equal

    % ------- Matrix Assembly -------
    % Extend element stiffness matrices to include rotational DOFs (phi_z unused)
    Kelmnow = [Kelm     zeros(size(Kelm,1), size(Kelm,1)/5);
               zeros(size(Kelm,1)/5, size(Kelm,1)*6/5)];

    % Rotate and assemble into global stiffness matrix
    Ka(indelm, indelm) = Ka(indelm, indelm) + Tnow2 \ Kelmnow * Tnow2;

    % Extend element mass matrices to include rotational DOFs
    Melmnow = [Melm     zeros(size(Kelm,1), size(Kelm,1)/5);
               zeros(size(Kelm,1)/5, size(Kelm,1)*6/5)];

    % Rotate and assemble into global mass matrix
    Ma(indelm, indelm) = Ma(indelm, indelm) + Tnow2 \ Melmnow * Tnow2;

end

% Identify and remove unused DOFs (e.g., phi_z DOFs not present in model)
indF = (find(sum(abs(Ma))==0))';
Ka(indF,:) = [];
Ka(:,indF) = [];
Ma(indF,:) = [];
Ma(:,indF) = [];
sem_mesh.ind_ALL(indF,:) = [];
sem_mesh.posn(indF,:) = [];
clear indF

% Report assembly time
disp(['Assembly: ' num2str(round(toc,1)) ' s'])

%% ------------------------------------------------------------------------
% -------------- Boundary Condition Application ---------------------------
% -------------------------------------------------------------------------
% This section applies the prescribed boundary conditions (BC) to the 
% assembled global stiffness (Ka) and mass (Ma) matrices. The function 
% modifies the system by:
%   - Identifying degrees of freedom (DOFs) that are constrained (e.g., fixed, simply supported)
%   - Applying coordinate transformations where necessary (e.g., for rotational DOFs)
%   - Eliminating constrained DOFs from the global system
%   - Updating the mesh structure (sem_mesh) to reflect removed DOFs
% The resulting system matrices (Ka_mod, Ma_mod) are ready for solving
% free vibration or forced response problems with the imposed BCs.

tic

[Ka_mod, Ma_mod, sem_mesh] = Boundary_Conditions_Application(Ka, Ma,...
    sem_mesh, BC);

% Report boundary application time
disp(['Boundary Condition: ' num2str(round(toc,1)) ' s'])

%% ------------------------------------------------------------------------
% -------------- Eigenvalue Solution --------------------------------------
% -------------------------------------------------------------------------
% This section solves the generalized eigenvalue problem for the modified
% system matrices (after boundary condition application). The goal is to 
% compute the first few natural frequencies and their corresponding mode shapes.

tic

% Shift value (sigma) for improved numerical stability and convergence
% Used with the shift-invert mode in `eigs` to target low-frequency modes
sigma = 0.01; 

% Compute the 20 smallest eigenvalues and eigenvectors of the system:
%   Ka_mod * x = lambda * Ma_mod * x
% The 'eigs' function uses sparse matrix methods for efficiency.
[eigVec,eigVal,flag] = eigs(Ka_mod,Ma_mod,20,sigma);

% Extract natural frequencies (rad/s) by taking square roots of eigenvalues
% The sigma shift is subtracted back from the eigenvalues before square rooting
[wns, loc] = sort(real(sqrt(diag(eigVal) - sigma)));

% Reorder eigenvectors accordingly
eigVec = eigVec(:, loc);

% Display the first 20 natural frequencies in Hz
disp((wns(1:20))/2/pi)

% Display elapsed time for eigenvalue solution
disp(['Eigen Solution time: ' num2str(round(toc,1)) ' s'])