function [A, B, D, G] = Stacking_LP_ABD(theta_vec, matinv, h)

    theta = theta_vec*pi/180;   % Fiber angles in degree
    NN = length(theta);
    
    % Lamination parameters for in-plane stiffness matrix
    
    V1A = 0;
    V2A = 0;
    V3A = 0;
    V4A = 0;
    
    for i = 1:1:NN
        
        V1A = (1/NN)*(cos(2*theta(i)))+V1A;
        V2A = (1/NN)*(sin(2*theta(i)))+V2A;
        V3A = (1/NN)*(cos(4*theta(i)))+V3A;
        V4A = (1/NN)*(sin(4*theta(i)))+V4A;
        
    end
    
    % Coupling lamination parameters
    
    
    V1B = 0;
    V2B = 0;
    V3B = 0;
    V4B = 0;
    
    for i= 1:1:NN
        
        V1B = (2/NN^2)*(cos(2*theta(i)))*((-NN/2+i)^2-(-NN/2+i-1)^2)+V1B;
        V2B = (2/NN^2)*(sin(2*theta(i)))*((-NN/2+i)^2-(-NN/2+i-1)^2)+V2B;
        V3B = (2/NN^2)*(cos(4*theta(i)))*((-NN/2+i)^2-(-NN/2+i-1)^2)+V3B;
        V4B = (2/NN^2)*(sin(4*theta(i)))*((-NN/2+i)^2-(-NN/2+i-1)^2)+V4B;
        
    end
    
    %%Bending lamination parameters
    
    V1D = 0;
    V2D = 0;
    V3D = 0;
    V4D = 0;
    
    for i= 1:1:NN
        
        V1D = (4/NN^3)*(cos(2*theta(i)))*((-NN/2+i)^3-(-NN/2+i-1)^3)+V1D;
        V2D = (4/NN^3)*(sin(2*theta(i)))*((-NN/2+i)^3-(-NN/2+i-1)^3)+V2D;
        V3D = (4/NN^3)*(cos(4*theta(i)))*((-NN/2+i)^3-(-NN/2+i-1)^3)+V3D;
        V4D = (4/NN^3)*(sin(4*theta(i)))*((-NN/2+i)^3-(-NN/2+i-1)^3)+V4D;
        
    end
    
    
    % In-Plane Stiffness matrix
    
    [VA] = [1 V1A V3A 0 0 ;
        1 -V1A V3A 0 0 ;
        0 0 -V3A 1 0 ;
        0 0 -V3A 0 1 ;
        0 V2A/2 V4A 0 0;
        0 V2A/2 -V4A 0 0];
    
    %% Shear stiffness matrices
    
    [GA] = [ 1 -V1A;
        1 V1A;
        0 -V2A];
    
    U_invariant = [matinv.U1 matinv.U2 matinv.U3 matinv.U4 matinv.U5]';
    U_shear = [matinv.U11 matinv.U22]';
    
    %% In-plane stiffness matirices
    
    A = h*VA*U_invariant;                   %% In-plane laminate stiffness for core
    
    
    %% Transverse Shear stiffness matrices
    
    G = h*GA*U_shear;                  %% Transverse shear stiffness matrices for core
    
    % Coupling Stiffness matrix
    
    [VB] = [0 V1B V3B 0 0;
        0 -V1B V3B 0 0;
        0 0 -V3B 0 0;
        0 0 -V3B 0 0;
        0 V2B/2 V4B 0 0;
        0 V2B/2 -V4B 0 0];
    
    
    B = h^2/4*VB*U_invariant;
    
    
    % Bending Stiffness Matrix
    
    [VDD] = [1 V1D V3D 0 0;
        1 -V1D V3D 0 0;
        0 0 -V3D 1 0;
        0 0 -V3D 0 1;
        0 V2D/2 V4D 0 0;
        0 V2D/2 -V4D 0 0];
    
    
    D = h^3/12*VDD*U_invariant;
    

end