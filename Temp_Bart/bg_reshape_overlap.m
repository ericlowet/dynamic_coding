function [Xres, winStart]=bg_reshape_overlap(X,dWin,winLen,dim)
% Xres=bg_reshape_overlap(X,dWin,winLen,dim)
%
% Reshapes matrix along dimension "dim", but with new overlapping  windows
% defined by dWin and winLen;
% 
% useful for Welch's method FT

if nargin<4
  dim=1;
end

tLen=size(X,dim);
winStart=1:dWin:(tLen-winLen+1);
numWin=numel(winStart);


nd=ndims(X);

%move time dimension to first dimension
if dim~=1
  perm=[1:nd];
  sel=true(size(perm));
  sel(dim)=false;
  perm=[dim perm(sel)];
  X=permute(X,perm);
end

szX=size(X);



% collapse all non-time dimensions;
X=reshape(X,tLen,[]);
Xres=nan(winLen,numel(winStart),prod(szX(2:end)));
winSel=repmat([1:winLen]',1,numel(winStart));
winSel=bsxfun(@plus,winSel,winStart-1);
for n=1:size(X,2)
  Xres(:,:,n)=reshape(X(winSel(:),n),size(winSel));
end

% reshape back all the non-time dimenions
Xres=reshape(Xres, [winLen, numWin, szX(2:end)]);

% move window number dimension to back
Xres=permute(Xres, [1 3:numel(size(Xres)) 2]);

if dim~=1
  % go back to original dimension order (with added last dimension)
  Xres=ipermute(Xres, [perm nd+1]);
end