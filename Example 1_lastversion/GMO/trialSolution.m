function [y, dy, d2y] = trialSolution(model,x,a,b,A,B)
    %n = length(x)-1;
    n = length(x);
    h = 0.01;%(b-a)/n;
    Xa = x-a;
    Xb = x-b;
    XaXb = Xa.*Xb;
    X = [x x+h x-h];
    Y = model(X);
    %net = Y(1:n);
    %dnet = (model(x+h) + model(x-h))/(2*h);
    dnet = (Y(n+1:2*n) + Y(2*n+1:3*n))/(2*h);
    %d2net = (model(x+h) - 2*net + model(x-h))/(h^2);
    d2net = (Y(n+1:2*n)- 2*Y(1:n) +  Y(2*n+1:3*n))/(h^2);
    
    %     y = B*(x - a)./ (b-a) - A*(x-b)./(b-a) + (x-a).*(x-b).*net;
    %     dy = (B-A)/(b-a) + (2*x-a-b).*net +(x-a).*(x-b).*dnet;
    %     d2y = 2*net + 2*(2*x-a-b).*dnet +(x-a).*(x-b).*d2net;
    y = B*Xa./ (b-a) - A*Xb./(b-a) + XaXb.*Y(1:n);
    dy = (B-A)/(b-a) + (Xa + Xb).*Y(1:n) + XaXb.*dnet;
    d2y = 2*Y(1:n) + 2*(Xa + Xb).*dnet + XaXb.*d2net;
end