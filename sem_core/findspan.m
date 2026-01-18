function span = findspan(n,p,u,U)

%The NURBS Book, p.68, Algorithm A2.1

if u == U(n+2)
    span = n;
    return;
end;

low = p;
high = n+1;
mid = floor((low+high)/2);

while u<U(mid+1) || u>=U(mid+2);
    if u < U(mid+1);
        high = mid;
    else
        low = mid;
    end;
    mid = floor((low+high)/2);
end;

span = mid;
span = span + 1;