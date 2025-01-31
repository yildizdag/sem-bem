function [x] = lobat(n)

% -------------------------------------------------------------------------
% This function calculates G/L points spanning -1 to 1, given n the number 
% of sampling points required
% -------------------------------------------------------------------------

nm1=n-1 ;
x=sin(pi*(-nm1+2*(0:nm1))/(2*nm1));
x=x';
