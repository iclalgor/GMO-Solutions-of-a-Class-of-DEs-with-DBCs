function y = N(x,alpha,w,bias)
 n = length(x);
 m =length(alpha);
 y = zeros(1,n);
 z = zeros(n,m);
 for j=1:n
   for i = 1:m
     z(j,i) = w(i)*x(j)+bias(i);
     y(j) = y(j) + alpha(i)*sigma(z(j,i));
   end
 end
end