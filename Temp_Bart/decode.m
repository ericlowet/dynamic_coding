function z=decode(rate,theta,X,Y,numOrient)
% z=decode(rate,theta,X,Y,numOrient)

% x=unique(X);
x=min(X(:)):max(X(:));
% y=unique(Y);
y=min(Y(:)):max(Y(:));
theta=theta/numOrient;
theta=theta*2*pi;
z=nan(max(y),max(x),size(rate,2));
% z=nan(numel(y),numel(x));

nansel=isnan(rate);
rate(nansel)=0;
for nx=1:numel(x)
  for ny=1:numel(y)    
    selXY=X==x(nx) & Y==y(ny);
    
    z(ny,nx,:)=(rate(selXY,:)'*exp(1i*theta(selXY)))./(sum(rate(selXY,:)).');
    z(ny,nx,any(nansel(selXY,:)))=nan;
  end
end 
