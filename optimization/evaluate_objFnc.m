function obj = evaluate_objFnc(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp)
%--------------------------------------------------------------------------
% Update Geometry
[Nurbs2D_plate,Nurbs2D_stiff] = update_baseline(Nurbs2D_plate,Nurbs2D_stiff,pconn,dcp);
%--------------------------------------------------------------------------
figure;
hold on
iga2DmeshPlotNURBS(Nurbs2D_plate);
iga2DmeshPlotNURBS(Nurbs2D_stiff);
hold off
%-Young's Modulus
E = 205E9;
nu = 0.3;
rho = 7800;
%-Geometric Props
t = 0.004;   %-thickness
%-Number of Tchebychev Polynomials (per element)
N = 5;
%-Element Type:
ET = 2; % 1: Plate on x-y plane (3 DOF)
        % 2: Shell in 3D (6 DOF)
%-Formulation:
form = 2; % 1: Based on NURBS
          % 2: Based on Chebyshev
%-DOF per Sampling Point:
if ET == 1
    shell_dof = 3;
elseif ET == 2
    shell_dof = 6;
end
% Create SEM Mesh
semOpt2D = semOpt2Dmesh(Nurbs2D_plate,Nurbs2D_stiff,N,shell_dof);
%--------------------------------------------------------------------------
semOpt2D.ET = ET;
semOpt2D.form = form;
semOpt2D.N = N;
semOpt2D.non = size(semOpt2D.nodes,1);
semOpt2D.dof = semOpt2D.non*shell_dof;
semOpt2D.E = E;
semOpt2D.nu = nu;
semOpt2D.rho = rho;
semOpt2D.t = t;
semOpt2D.D = (E*t^3)/12/(1-nu^2);
semOpt2D.G = E/2/(1+nu);
semOpt2D.Ds = (5/6)*semOpt2D.G*t;
semOpt2D.lame = 2*semOpt2D.G/(1-semOpt2D.nu);  
%
[K,M] = global2D(semOpt2D);
%
x_min = min(semOpt2D.nodes(:,1)); x_max = max(semOpt2D.nodes(:,1));
y_min = min(semOpt2D.nodes(:,2)); y_max = max(semOpt2D.nodes(:,2));
ind = find(semOpt2D.nodes(:,1)<x_min+1E-6 | semOpt2D.nodes(:,1)>x_max-1E-6 |...
           semOpt2D.nodes(:,2)<y_min+1E-6 | semOpt2D.nodes(:,2)>y_max-1E-6);
BounNodes = unique([6.*ind-5; 6.*ind-4; 6.*ind-3; 6.*ind-2; 6.*ind-1; 6.*ind]);
%
K(BounNodes,:) = []; K(:,BounNodes) = [];
M(BounNodes,:) = []; M(:,BounNodes) = [];
%
tic;
%-Eigenvalue Solver
sigma = 0.1;
[V,freq] = eigs(K,M,1,sigma);
[freq,loc] = sort((sqrt(diag(freq)-sigma)));
%
%V = V(:,loc);
obj = freq/2/pi;
end