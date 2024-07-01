function E = Cost(model,x,a,b,A,B,W)
  model = setwb(model,W);
  [y,dy,d2y] = trialSolution(model,x,a,b,A,B);
  F = f(x,y,dy);
  E = sum((d2y - F).^2)/length(x);
end
