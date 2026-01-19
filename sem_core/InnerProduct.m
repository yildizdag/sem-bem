function V = InnerProduct(space)


    %% Order is increased due to the product
    
    % The order of the obtained function is the total 
    % of the functions producted.
    space2 = space;
    space2.N = 2*space.N; 
    
    [FT2,IT2] = cheb(space2);
    [FT,~] = cheb(space); 
    
    %% Calculation of the Inner Product Matrix
    
    % The Tchebychev definite integral matrix
    
    v_di = cheb_di(space2);
    v_d2N = diag(v_di'*FT2); % Appendix D.2 in Yagci et al. 
    
    I = diag(ones(space.N,1));
    Z = 0*I;
    S = IT2*[I;Z]*FT;
    
    V = S'*v_d2N*S;

end