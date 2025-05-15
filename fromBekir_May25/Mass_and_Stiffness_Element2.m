function [xelm,yelm,Kelm,Melm] = ...
    Mass_and_Stiffness_Element2(ro, E, pois, Lxi, Leta, h, polynums, loc)

    polynum_xi = polynums(1);
    polynum_eta = polynums(2);

    % Discretization in each direction and calculation of parameters in mapped frame
    [~,~,~,xi,V_xi,Q1_xi,~,space_xi] = Discretization(Lxi, polynum_xi,'xi');
    [~,~,~,eta,V_eta,Q1_eta,~,space_eta] = Discretization(Leta, polynum_eta,'eta');
    
    Mapping_Order = 4;
    [xelm, yelm, dxdxi, dydxi, dxdeta, dydeta, ~, ~] = ...
        Cross_section_Mapping(Mapping_Order, loc, xi, eta);
    xelm = xelm';
    yelm = yelm';
    
    % Checking the Jacobian and derivative ralations between the physical frame and mapped frame:
    %       - D: Jacobian for the mapping
    dxidx = zeros(length(xi),length(eta));
    detadx = zeros(length(xi),length(eta));
    dxidy = zeros(length(xi),length(eta));
    detady = zeros(length(xi),length(eta));
    JAC = zeros(length(xi),length(eta));
    
    for i=1:length(xi)
        for j=1:length(eta)
            Dij = [dxdxi(i,j)   dydxi(i,j);...
                   dxdeta(i,j)  dydeta(i,j)];
            JAC(i,j) = det(Dij);
            Eij = inv(Dij);
    
            dxidx(i,j) = Eij(1,1);
            detadx(i,j) = Eij(1,2);
    
            dxidy(i,j) = Eij(2,1);
            detady(i,j) = Eij(2,2);
    
        end
    end
    
    %min(min(JAC))
    
    % Applying the vector mapping
    [VD, QDx, QDy] = Vector_mapping_2D_general (polynum_xi, polynum_eta, ...
                  xi, eta, V_xi, V_eta, Q1_xi, Q1_eta, JAC, ...
                  dxidx, detadx, dxidy, detady, space_xi, space_eta);
     
    mu = E/2/(1+pois);              % Shear Modulus
    lame = 2*mu/(1-pois);           % Lame Constant
    k_sc = sqrt(pi^2/12);           % sqrt(pi^2/12);
    beta_1 = 0;                     % curvature rate
    beta_2 = 0;                     % curvature rate
    
    % B) Calculation of system matrices
    [Kelm, Melm] = Mass_and_Stiffness_FSDT_DC(polynum_xi, polynum_eta, ...
        ro, lame, mu, pois, k_sc, h,...
        VD, QDx, QDy, beta_1, beta_2);



end