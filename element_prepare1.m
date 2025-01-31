function [locs,xlocalnow,ylocalnow] = ...
    element_prepare1(elementnow,nodes)

xlocalnow = [1 0 0];
ylocalnow = [0 1 0];

element_renum = elementnow([1 5 25 21 2 10 24 16 3 15 23 11 4 20 22 6 7 9 19 17 8 14 18 12 13]);

locs = [nodes(element_renum,1) nodes(element_renum,2)];




