function Nurbs2D = iga2Dmesh(fileName,numpatch,local_dof)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS FUNCTION CREATES IGA DATA STRUCTURE 
% WITH THE GIVEN INPUT FILES FOR 2D PROBLEMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileName: Name of the input file group
% numpatch: Number of Patch to Read
% local_dof: Number of DOF per each control point (node)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nurbs2D.numpatch = numpatch;
Nurbs2D.local_dof = local_dof;
count = 1;
for i = 1:numpatch
    file = [fileName num2str(count)];
    fileID = fopen(file);
    patchInfo = fscanf(fileID,'%f');
    %Get Knot vectors:
    vKnotCount = patchInfo(1);
    Vknots = patchInfo(2:(vKnotCount+1));
    uKnotCount = patchInfo(vKnotCount+2);
    Uknots = patchInfo(vKnotCount+3:(vKnotCount+uKnotCount+2));
    Nurbs2D.knots.V{i} = Vknots;
    Nurbs2D.knots.U{i} = Uknots;
    %Get Control Points Net
    pointCount = patchInfo(vKnotCount+uKnotCount+3);
    pointStart = vKnotCount+uKnotCount+4;
    cPoints = reshape(patchInfo(pointStart:end),4,pointCount);
    Nurbs2D.nodes{i} = transpose(cPoints);
    %Order:
    [~, ia1, ~] = unique(Uknots,'last');
    [~, ia2, ~] = unique(Vknots,'last');
    Uorder = ia1(1); Vorder = ia2(1);
    Nurbs2D.order{i}(1,1) = Uorder; Nurbs2D.order{i}(2,1) = Vorder;
    %Number of Basis Function:
    Unumber = numel(Uknots)-Uorder; Vnumber = numel(Vknots)-Vorder;
    Nurbs2D.number{i}(1,1) = Unumber; Nurbs2D.number{i}(1,2) = Vnumber;
    %Control Points:
    cPoints = reshape(cPoints,4,Unumber,Vnumber);
    Nurbs2D.cPoints{i} = cPoints;
    %Get Patch Boundary Data:
    Nurbs2D.edges.knots.U{i,1} = Uknots; Nurbs2D.edges.knots.U{i,2} = Vknots;
    Nurbs2D.edges.knots.U{i,3} = Uknots; Nurbs2D.edges.knots.U{i,4} = Vknots;
    Nurbs2D.edges.coefsNo{i,1} = 1:Unumber; Nurbs2D.edges.coefsNo{i,2} = Unumber:Unumber:(Vnumber*Unumber); 
    Nurbs2D.edges.coefsNo{i,3} = (Vnumber*Unumber):-1:(Vnumber*Unumber-Unumber+1); Nurbs2D.edges.coefsNo{i,4} = (Vnumber*Unumber-Unumber+1):-Unumber:1;
    Nurbs2D.edges.order{i,1} = Uorder; Nurbs2D.edges.order{i,2} = Vorder;
    Nurbs2D.edges.order{i,3} = Uorder; Nurbs2D.edges.order{i,4} = Vorder;
    Nurbs2D.edges.number{i,1} = Unumber; Nurbs2D.edges.number{i,2} = Vnumber;
    Nurbs2D.edges.number{i,3} = Unumber; Nurbs2D.edges.number{i,4} = Vnumber;
    count = count + 1;
end

for i = 1:Nurbs2D.numpatch
    
    [INC,IEN,nel,nnp,nen] = connectivity(Nurbs2D.order{i},Nurbs2D.number{i});
    Nurbs2D.INC{i} = INC;
    Nurbs2D.IEN{i} = IEN;
    Nurbs2D.nel{i} = nel;
    Nurbs2D.nnp{i} = nnp;
    Nurbs2D.nen{i} = nen;
    
    for j = 1:4
        [INC,IEN,nel,nnp,nen] = connectivity(Nurbs2D.edges.order{i,j},Nurbs2D.edges.number{i,j});
        Nurbs2D.edges.INC{i,j} = INC;
        Nurbs2D.edges.IEN{i,j} = IEN;
        Nurbs2D.edges.nel{i,j} = nel;
        Nurbs2D.edges.nnp{i,j} = nnp;
        Nurbs2D.edges.nen{i,j} = nen;
    end
    
    el1=Nurbs2D.edges.number{i,1} - Nurbs2D.edges.order{i,1} + 1;
    el2=Nurbs2D.edges.number{i,2} - Nurbs2D.edges.order{i,2} + 1;
    Nurbs2D.edges.el{i,1} = 1:el1;
    Nurbs2D.edges.el{i,2} = el1 : el1 : el1*el2;
    Nurbs2D.edges.el{i,3} = el1*el2 : -1 : el1*(el2-1)+1;
    Nurbs2D.edges.el{i,4} = el1*(el2-1)+1 : -el1 : 1;
    
    %Nurbs2D.int{i} = gauss2d(3,3);
    %[xgp,wgp,ngp] = gaussQuad2d(6,6);
    %Nurbs2D.xgp{i} = xgp;
    %Nurbs2D.wgp{i} = wgp;
    %Nurbs2D.ngp{i} = ngp;
    
end

%ID
eqn = 0; %Start counting degrees of freedoms
for i = 1:Nurbs2D.numpatch
    NNP = Nurbs2D.nnp{i};
    ID = zeros(Nurbs2D.local_dof,NNP);
    for j = 1:NNP
        for k = 1:Nurbs2D.local_dof
            eqn = eqn+1;
            ID(k,j) = eqn;
        end
    end
    Nurbs2D.ID{i} = ID;
end
Nurbs2D.eqn = eqn;

%LM
for i = 1:Nurbs2D.numpatch
    LM = zeros(Nurbs2D.local_dof,Nurbs2D.nen{i},Nurbs2D.nel{i});
    for j = 1:Nurbs2D.nel{i}
        for k = 1:Nurbs2D.nen{i}
            LM(:,k,j) = Nurbs2D.ID{i}(:, Nurbs2D.IEN{i}(k,j));
        end
    end
    Nurbs2D.LM{i} = LM;
end

%Boundary Connectivity
for i = 1:Nurbs2D.numpatch
    for j = 1:4
        Nurbs2D.edges.ID{i,j} = Nurbs2D.ID{i}(1:Nurbs2D.local_dof,Nurbs2D.edges.coefsNo{i,j});
        Nurbs2D.edges.nodes{i,j} = Nurbs2D.nodes{i}(Nurbs2D.edges.coefsNo{i,j},:);
    end
end

for i = 1:Nurbs2D.numpatch
    for j = 1:4
        Nurbs2D.edges.LM{i,j}=zeros(Nurbs2D.local_dof,Nurbs2D.edges.nen{i,j},Nurbs2D.edges.nel{i,j});
        for k = 1:Nurbs2D.edges.nel{i,j}
            for l = 1:Nurbs2D.edges.nen{i,j}
                Nurbs2D.edges.LM {i,j}(:,l,k) = Nurbs2D.edges.ID{i,j}(:, Nurbs2D.edges.IEN{i,j}(l,k));
            end
        end
    end
end