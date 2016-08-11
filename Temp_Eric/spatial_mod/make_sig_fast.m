function  [ signals, selm ] = make_sig(rastersp,R_e,net_size)


%[ signals ] = make_sig(rastersp(1:6400,:),R_e,net_size)
 
 
gridSize=net_size.*50e-1;
pdfConEE=single(connec_prob(0.6*50e-1,gridSize,R_e)>0.0001);

selm=find(pdfConEE(780,:)==1);

 t=zeros(net_size);
t(1:1:end,1:1:end)=1; %figure,imagesc(t)
selvect=find(t(:)==1); %length(selvect)
%rastersp=rastersp';
% figure,imagesc(pdfConEE)
rastersp=single(rastersp(:,1:2:end));
signals=zeros(length(selvect),size(rastersp,2),'single');
rastersp=rastersp';
timla=size(rastersp,1);
wind2=500;
tr=floor(timla/wind2);

n=0;      fs=2000;
for ind22=selvect'
    n=n+1
y=[];
for trl=1:tr
  [y1] =((rastersp((1:wind2)+wind2*(trl-1),:)*pdfConEE(:,ind22)));%, fs, [19 60]); 
  y =[y ; y1 ];
end
 [y1] =((rastersp(wind2*(trl)+1:end,:)*pdfConEE(:,ind22)));%, fs, [19 60]); 
  y =[y ; y1 ];

signals(n,:)=((y) );


end



%figure,imagesc(angle(signals))





function [pdfCon]=connec_prob(sig, gridSize, R_send, R_receive)

if nargin<4
  R_receive=R_send;
  selfFlag=0;
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
