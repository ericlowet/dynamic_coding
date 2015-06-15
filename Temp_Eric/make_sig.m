function  [ signals ] = make_sig(rastersp,R_e,net_size)


%[ signals ] = make_sig(rastersp(1:6400,:),R_e,net_size)
 
 
gridSize=net_size.*50e-6;
pdfConEE=connec_prob(1*50e-6,gridSize,R_e);

 t=zeros(net_size);
t(1:1:end,1:1:end)=1; %figure,imagesc(t)
selvect=find(t(:)==1); %length(selvect)

% figure,imagesc(pdfConEE)
signals=zeros(length(selvect),size(rastersp,2),'single');
n=0;      fs=2000;
for ind22=selvect'

    n=n+1;
       fs=2000;
[y] = ft_preproc_bandpassfilter(fastsmooth(rastersp'*pdfConEE(:,ind22),3,1,1)', fs, [22 70]);
%  [y] =(fastsmooth(rastersp'*pdfConEE(:,ind22),3,1,1)');%, fs, [19 60]);
  
signals(n,:)=hilbert(y);


end



%figure,imagesc(angle(signals))





function [pdfCon]=connec_prob(sig, gridSize, R_send, R_receive)

if nargin<4
  R_receive=R_send;
  selfFlag=1;
else
  selfFlag=0;
end

x_send=R_send(:,1);
y_send=R_send(:,2);

x_receive=R_receive(:,1);
y_receive=R_receive(:,2);

DX=abs(bsxfun(@minus,x_receive,x_send'));
DY=abs(bsxfun(@minus,y_receive,y_send'));


if numel(gridSize)==1
  gridSize=[gridSize gridSize];
end
DX=min(DX,gridSize(1)-DX);
DY=min(DY,gridSize(2)-DY);


% gaussian connection pdf
pdfX=exp(-DX.^2/(2*sig^2));
pdfY=exp(-DY.^2/(2*sig^2));
pdfCon=pdfX.*pdfY;

if selfFlag
  pdfCon(logical(eye(size(pdfCon))))=0; % no self connections
end
pdfCon=bsxfun(@rdivide, pdfCon, sum(pdfCon,2));
end

end
