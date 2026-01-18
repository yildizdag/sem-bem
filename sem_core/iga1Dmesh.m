function Nurbs1D = iga1Dmesh(fileName,numpatch,local_dof)
% Create 1-D NURBS Data for IGA
% Suitable for both Single- or Multi-patch Representations
Nurbs1D.numpatch = numpatch;
Nurbs1D.local_dof = local_dof;
count = 1;
for i = 1:numpatch
    file = [fileName num2str(count)];
    fileID = fopen(file);
    patchInfo = fscanf(fileID,'%f');
    %Get Knot vector:
    uKnotCount = patchInfo(1);
    Uknots = patchInfo(2:(uKnotCount+1));
    Nurbs1D.knots.U{i} = Uknots;
    %Get Control Points Net:
    pointCount = patchInfo(uKnotCount+2);
    pointStart = uKnotCount+3;
    cPoints = reshape(patchInfo(pointStart:end),4,pointCount);
    Nurbs1D.cPoints{i} = cPoints;
    %Order:
    [~, ia1, ~] = unique(Uknots,'last');
    Uorder = ia1(1);
    Nurbs1D.order{i}(1,1) = Uorder;
    %Number of Basis Function:
    Unumber = numel(Uknots)-Uorder;
    Nurbs1D.number{i}(1,1) = Unumber;
    %
    count = count + 1;
end
%
for i = 1:Nurbs1D.numpatch
    %
    [INC,IEN,nel,nnp,nen] = connectivity(Nurbs1D.order{i},Nurbs1D.number{i});
    Nurbs1D.INC{i} = INC;
    Nurbs1D.IEN{i} = IEN;
    Nurbs1D.nel{i} = nel;
    Nurbs1D.nnp{i} = nnp;
    Nurbs1D.nen{i} = nen;
end
%ID
eqn = 0;
for i = 1:Nurbs1D.numpatch
    NNP = Nurbs1D.nnp{i};
    ID = zeros(Nurbs1D.local_dof,NNP);
    for j = 1:NNP
        for k = 1:Nurbs1D.local_dof
            eqn = eqn+1;
            ID(k,j) = eqn;
        end
    end
    Nurbs1D.ID{i} = ID;
end
Nurbs1D.eqn = eqn;
%LM
for i = 1:Nurbs1D.numpatch
    LM = zeros(Nurbs1D.local_dof,Nurbs1D.nen{i},Nurbs1D.nel{i});
    for j = 1:Nurbs1D.nel{i}
        for k = 1:Nurbs1D.nen{i}
            LM(:,k,j) = Nurbs1D.ID{i}(:, Nurbs1D.IEN{i}(k,j));
        end
    end
    Nurbs1D.LM{i} = LM;
end