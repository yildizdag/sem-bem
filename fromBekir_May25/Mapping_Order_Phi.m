function [psixi, psieta, constants, locA] = Mapping_Order_Phi (Mapping_Order)

%------------------------------------------------------------------------
%
%
%------------------------------------------------------------------------

switch Mapping_Order
    case 1
        disp ('daha girilmedi :D');
    case 2
        
        locA=1:1:9; 
        
        psixi(1,:)=conv([1 0],[1 -1]);
        psixi(2,:)=conv([1 0],[1 1]);
        psixi(3,:)=conv([1 0],[1 1]);  
        psixi(4,:)=conv([1 0],[1 -1]);  
        
        psixi(5,:)=conv([1 -1],[1 1]);
        psixi(6,:)=conv([1 0],[1 1]);        
        psixi(7,:)=conv([1 -1],[1 1]);
        psixi(8,:)=conv([1 0],[1 -1]);
        
        psixi(9,:)=conv([1 -1],[1 1]);
        
        psieta(1,:)=conv([1 0],[1 -1]);
        psieta(2,:)=conv([1 0],[1 -1]);
        psieta(3,:)=conv([1 0],[1 1]);  
        psieta(4,:)=conv([1 0],[1 1]);     
        
        psieta(5,:)=conv([1 0],[1 -1]);
        psieta(6,:)=conv([1 -1],[1 1]);        
        psieta(7,:)=conv([1 0],[1 1]);
        psieta(8,:)=conv([1 -1],[1 1]);
        
        psieta(9,:)=conv([1 -1],[1 1]);
        
        constants=[1/4*ones(4,1); -1/2*ones(4,1); 1]; 
        
    case 3
        disp ('daha girilmedi :D');        
    case 4
        
        locA=1:1:25; 
        
        psixi(1,:)=conv([1 -1 0],[1 0 -1/4]);
        psixi(2,:)=conv([1 1 0],[1 0 -1/4]);
        psixi(3,:)=conv([1 1 0],[1 0 -1/4]);
        psixi(4,:)=conv([1 -1 0],[1 0 -1/4]);
        
        psixi(5,:)=conv([1 -1/2 0],[1 0 -1]);
        psixi(6,:)=conv([1 1 0],[1 0 -1/4]);
        psixi(7,:)=conv([1 1/2 0],[1 0 -1]);
        psixi(8,:)=conv([1 -1 0],[1 0 -1/4]);

        psixi(9,:)=conv([1 0 -1],[1 0 -1/4]);
        psixi(10,:)=conv([1 1 0],[1 0 -1/4]);
        psixi(11,:)=conv([1 0 -1],[1 0 -1/4]);
        psixi(12,:)=conv([1 -1 0],[1 0 -1/4]);

        psixi(13,:)=conv([1 1/2 0],[1 0 -1]);
        psixi(14,:)=conv([1 1 0],[1 0 -1/4]);
        psixi(15,:)=conv([1 -1/2 0],[1 0 -1]);
        psixi(16,:)=conv([1 -1 0],[1 0 -1/4]);

        psixi(17,:)=conv([1 -1/2 0],[1 0 -1]);
        psixi(18,:)=conv([1 1/2 0],[1 0 -1]);
        psixi(19,:)=conv([1 1/2 0],[1 0 -1]);
        psixi(20,:)=conv([1 -1/2 0],[1 0 -1]);

        psixi(21,:)=conv([1 0 -1],[1 0 -1/4]);
        psixi(22,:)=conv([1 1/2 0],[1 0 -1]);
        psixi(23,:)=conv([1 0 -1],[1 0 -1/4]);
        psixi(24,:)=conv([1 -1/2 0],[1 0 -1]);

        psixi(25,:)=conv([1 0 -1],[1 0 -1/4]);  
        
        psieta(1,:)=conv([1 -1 0],[1 0 -1/4]);
        psieta(2,:)=conv([1 -1 0],[1 0 -1/4]);
        psieta(3,:)=conv([1 1 0],[1 0 -1/4]);
        psieta(4,:)=conv([1 1 0],[1 0 -1/4]);

        psieta(5,:)=conv([1 -1 0],[1 0 -1/4]);
        psieta(6,:)=conv([1 0 -1],[1 -1/2 0]);
        psieta(7,:)=conv([1 1 0],[1 0 -1/4]);
        psieta(8,:)=conv([1 0 -1],[1 1/2 0]);

        psieta(9,:)=conv([1 -1 0],[1 0 -1/4]);
        psieta(10,:)=conv([1 0 -1],[1 0 -1/4]);
        psieta(11,:)=conv([1 1 0],[1 0 -1/4]);
        psieta(12,:)=conv([1 0 -1],[1 0 -1/4]);

        psieta(13,:)=conv([1 -1 0],[1 0 -1/4]);
        psieta(14,:)=conv([1 0 -1],[1 1/2 0]);
        psieta(15,:)=conv([1 1 0],[1 0 -1/4]);
        psieta(16,:)=conv([1 0 -1],[1 -1/2 0]);

        psieta(17,:)=conv([1 0 -1],[1 -1/2 0]);
        psieta(18,:)=conv([1 0 -1],[1 -1/2 0]);
        psieta(19,:)=conv([1 0 -1],[1 1/2 0]);
        psieta(20,:)=conv([1 0 -1],[1 1/2 0]);

        psieta(21,:)=conv([1 0 -1],[1 -1/2 0]);
        psieta(22,:)=conv([1 0 -1],[1 0 -1/4]);
        psieta(23,:)=conv([1 0 -1],[1 1/2 0]);
        psieta(24,:)=conv([1 0 -1],[1 0 -1/4]);

        psieta(25,:)=conv([1 0 -1],[1 0 -1/4]);

        constants=[4/9*ones(4,1); -16/9*ones(4,1); 8/3*ones(4,1); -16/9*ones(4,1); 64/9*ones(4,1); -32/3*ones(4,1); 16]; 

end