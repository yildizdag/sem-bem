function [u]=Interpol2D(spacex, spacey, x, y, a)

% Given an array of points x that lies between 
% xleft and xright, and an array of Tchebychev coefficients
% a that is used for function defined on the interval (xleft,xright)
% this returns a vector y giving the value of the 
% Tchebychev series at the points x

%% ----------------- Preallocation of vectors ---------------------
[rrx,ccx] = size(x);
if rrx < ccx
    x = transpose(x);
end

[rry,ccy] = size(y);
if rry < ccy
    y = transpose(y);
end

%% ----------------- Converting every domain to -1< <1 ------------

Nx = spacex.N; 
xleft= spacex.a;  
xright = spacex.b ;

alphax = (xright-xleft)/2;  
betax= (xright+xleft)/2 ;

zx = (x-betax)/alphax ;

Ny = spacey.N; 
yleft= spacey.a;  
yright = spacey.b ;

alphay = (yright-yleft)/2;  
betay= (yright+yleft)/2 ;

zy = (y-betay)/alphay ;


%% ----------------- Generating Tcheb Polys ---------------------

TMatx= TchebPoly(Nx-1,zx);
TMaty= TchebPoly(Ny-1,zy);

spnxy=Nx*Ny;
spnxyInt=length(x)*length(y);

TMAT=zeros(spnxyInt,spnxy);
for i=1:length(x)
    for j=1:length(y)  
        for i2=1:Nx
            for j2=1:Ny
                count=(i-1)*length(y)+j;
                count2=(i2-1)*Ny+j2;
                TMAT(count,count2)=TMAT(count,count2)+TMatx(i,i2)'*TMaty(j,j2)';
            end
        end
    end
end

u=TMAT*a;



