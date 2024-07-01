function yt = trialy(x,a,ya,b,yb,cf)
n = length(x);
m = length(cf)/3;
alpha = cf(1:m);
w = cf(m+1:2*m);
bias = cf(2*m+1:3*m);
yt =zeros(1,n);
for j=1:n
    yt(j) = (ya*(x(j) - b) - yb*(x(j) - a))/(a - b) +(x(j)-a)*(x(j)-b)*N(x(j),alpha,w,bias);
end

end