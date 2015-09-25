function z=decode(rate,theta,X,Y,numOrient)
% z=decode(rate,theta,X,Y,numOrient)

x=unique(X);
y=unique(Y);
theta=theta/numOrient;
theta=theta*2*pi;
z=nan(max(y),max(x));
% z=nan(numel(y),numel(x));

for nx=1:numel(x)
  for ny=1:numel(y)    
    selXY=X==x(nx) & Y==y(ny);
    
    z(ny,nx)=rate(selXY)'*exp(1i*theta(selXY))/(sum(rate(selXY)));
  end
end
