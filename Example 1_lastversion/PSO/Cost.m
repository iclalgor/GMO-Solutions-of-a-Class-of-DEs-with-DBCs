function y = Cost(position,x,xa,ya,xb,yb)
n = length(x);
m = length(position)/3;
alpha = position(1:m);
w = position(m+1:2*m);
bias = position(2*m+1:3*m);
dN_dx=zeros(1,n);
d2N_dx2=zeros(1,n);
Y=zeros(1,n);
DY=zeros(1,n);
D2Y=zeros(1,n);
z=zeros(1,m);
totalE = 0;
for j=1:n
    if (j==1)
        Y(1) = ya; DY(1) = 0;
    else
        Y(j) = (ya*(x(j) - xb) - yb*(x(j) - xa))/(xa - xb) +(x(j) - xa)*(x(j) - xb)*N(x(j),alpha,w,bias);
        for i=1:m
            z(i) = w(i)*x(j) + bias(i);
            dN_dx(j) = dN_dx(j) + alpha(i)*w(i)*sigma(z(i))*(1-sigma(z(i)));
            d2N_dx2(j) = d2N_dx2(j) + alpha(i)*(w(i))^2*sigma(z(i))*(1-sigma(z(i)))*(1-2*sigma(z(i)));
        end
        DY(j) = (ya - yb)/(xa - xb) + (2*x(j) - xa - xb)*N(x(j),alpha,w,bias) + (x(j) - xa)*(x(j) - xb)*dN_dx(j);
        D2Y(j) = 2*(N(x(j),alpha,w,bias) + (2*x(j) - xa - xb)*dN_dx(j)) + (x(j) - xa)*(x(j) - xb)*d2N_dx2(j);
    end
    totalE = totalE + (D2Y(j) - f(x(j),Y(j),DY(j)))^2;
    %totalE = totalE + abs(DY(j) - f(x(j),Y(j)));
end
y = totalE/n;
return
end