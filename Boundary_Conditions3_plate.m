function [Ka,Ma,indA,indB,posn,indAf] = Boundary_Conditions3_plate(Ka,Ma,indA,indB,BCs,posn)

indAf = [];

if BCs(1) == 'C'
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0))];
end
if BCs(2) == 'C'
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0))];
end
if BCs(3) == 'C'
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0))];
end
if BCs(4) == 'C'
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0))];
end

if BCs(1) == 'S'
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA == 1))];
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA == 2))];
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA == 3))];
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA == 4))];
end
if BCs(2) == 'S'
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA == 1))];
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA == 2))];
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA == 3))];
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA == 4))];
end
if BCs(3) == 'S'
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA == 1))];
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA == 2))];
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA == 3))];
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA == 5))];
end
if BCs(4) == 'S'
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA == 1))];
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA == 2))];
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA == 3))];
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA == 5))];
end


if BCs(1) == 'P'
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 1))];
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 2))];
    indAf = [indAf; find((posn(:,1) == min(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 3))];
end
if BCs(2) == 'P'
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 1))];
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 2))];
    indAf = [indAf; find((posn(:,1) == max(posn(:,1))).*(posn(:,3) == 0).*(indA(:,2) == 3))];
end
if BCs(3) == 'P'
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 1))];
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 2))];
    indAf = [indAf; find((posn(:,2) == min(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 3))];
end
if BCs(4) == 'P'
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 1))];
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 2))];
    indAf = [indAf; find((posn(:,2) == max(posn(:,2))).*(posn(:,3) == 0).*(indA(:,2) == 3))];
end

Ka(indAf,:) = [];
Ka(:,indAf) = [];
Ma(indAf,:) = [];
Ma(:,indAf) = [];

indA(indAf,:) = [];
indB(indAf,:) = [];
posn(indAf,:) = [];

end
