function z = f(x,y,Dy)
  z = (-(2.*x)./(1+x.^2)).*Dy+y+(2./(1+x.^2))-log(1+x.^2);
end