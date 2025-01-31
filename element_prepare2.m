function [indelm,Tnow2] = ...
    element_prepare2(xlocalnow,ylocalnow,elementpoints,indR)

% xlocalnow, ylocalnow: unit vector in local x- and y-directions
% zlocalnow, unit vector in local z-direction
zlocalnow = cross(xlocalnow,ylocalnow);


elementpoints = elementpoints(elementpoints>0);
indelm = [elementpoints'; elementpoints'+(length(indR)/6); 
    elementpoints'+(length(indR)/6)*2; elementpoints'+(length(indR)/6)*3; 
    elementpoints'+(length(indR)/6)*4; elementpoints'+(length(indR)/6)*5];

% Transformation matrix for system matrices, for a single sampling point
% Tnow1 is used if the sampling point is shared by multiple elements
% Tnow1no is used otherwise
Tnow1 = [xlocalnow 0 0 0;...
    ylocalnow 0 0 0;...
    zlocalnow 0 0 0;...
    0 0 0 -ylocalnow;...
    0 0 0 xlocalnow;...
    0 0 0 zlocalnow;];
Tnow1no = eye(6,6);

indRnow = indR(indelm);

% Transformation matrix for system matrices, for the assembly
Tnow2 = sparse(length(indelm),length(indelm));
for di2 = 1:length(elementpoints)
    if zlocalnow(3) == 1
        for di3 = 1:6
            for di4 = 1:6
                Tnow2((di3-1)*end/6+di2,(di4-1)*end/6+di2) = Tnow1(di3,di4);
            end
        end
    else
        if indRnow(di2) == 1
            for di3 = 1:6
                for di4 = 1:6
                    Tnow2((di3-1)*end/6+di2,(di4-1)*end/6+di2) = Tnow1(di3,di4);
                end
            end
        end
        if indRnow(di2) == 0
            for di3 = 1:6
                for di4 = 1:6
                    Tnow2((di3-1)*end/6+di2,(di4-1)*end/6+di2) = Tnow1no(di3,di4);
                end
            end
        end
    end
end






