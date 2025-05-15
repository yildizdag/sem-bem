function [x] = slobat(space)

% ------------------------------------------------------------------------
% This function calculates the physical sampling points in the domain of
% (xleft, xright) that correspond to Gauss-Lobatto (G/L) points in (-1,1). 
% The mapping is defined as, 
%
%                       x = alpha*zeta + beta
%
% Example: If we have f = sin(x) from 0 to pi, and we want to perform the
% sampling for 10 points, we would define 
%
%           space.N = 10; space.a = 0; space.b = pi;
%
% Then we would calculate 
%
%                     [xs] = slobat(space);
%
% We could now evaluate the sin(x) at points xs to sample it at G/L points
% to find underbar(f).
%
%    t = 0:pi/100:pi;
%    f1 = sin(t);
%    space.N = 10; space.a = 0; space.b=pi;
%    [x1] = slobat(space);
%    plot(t,sin(t),'b')
%    hold on
%    plot(x1,sin(x1),'r+')
%--------------------------------------------------------------------------

NTch = space.N;  xleft= space.a ;  xright = space.b ;

lobatto = lobat(NTch) ; %G/L points in (-1,1)
alpha = (xright-xleft)/2 ; 
beta = (xright+xleft)/2 ;
z = ones(1,NTch); 
x = alpha*lobatto + beta*z' ;

