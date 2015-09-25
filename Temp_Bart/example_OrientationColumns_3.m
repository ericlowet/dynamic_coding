cfg=[];
cfg.numOrient=8;
cfg.gridSz=4;
numE=50;
cfg.numE=numE;
% [nEE nEI nIE nII]
cfg.numConIntra=[10 15 20 20];
cfg.numConInter=[50 00 0 0];
cfg.sig_dist=[2 2 .5 0];
cfg.sig_theta=[.05 .05 .05 2]*pi;
cfg.strConIntra=[1e-3 5e-2 5e-2 2e-2];
cfg.strConInter=[.5e-3 3e-2 5e-2 2e-2];


[ conMat, idx, cfg, conMatSep] = genOrientColumns(cfg);


figure(8)
subaxis(1,3,1)
imagesc(conMat>0)
title('complete S')
subaxis(1,3,2)
imagesc(conMatSep{1,1}>0)
title('S_{intra}')
subaxis(1,3,3)
imagesc(conMatSep{1,2}>0)
title('S_{inter}')
colormap gray
%%
X=idx.X;
Y=idx.Y;
EI=idx.EI;
theta=idx.theta;
numNeur=numel(X);
numOrient=cfg.numOrient;
gridSz=cfg.gridSz;

sacFreq=3;
stimLen=10;
simLen=5e3;


% tuning functions
figure(7)
thh=linspace(-pi/2,pi/2,1001);
sig=pi/6;
TF=(cos(thh*2)+1)/2;
plot(thh,TF)
title('Tuning function')
ylabel('input strength')
xlabel('stimulus orientation')
axis tight



phi=-ceil([0:stimLen:simLen-stimLen]/1e3*sacFreq)*pi/4+0*pi/2;
% phi=cos([0:stimLen:simLen-stimLen]/1e3*2*pi*sacFreq*.65)*pi/2;
% phi=cumsum(pi/8*randn(numel([0:stimLen:simLen-stimLen]),1));
% phi=cumsum(sign(rand(numel([0:stimLen:simLen-stimLen]),1)-.5)*pi/32);


I=zeros(numNeur,numel(phi));
I(EI,:)=1;

tbins=0:stimLen:stimLen*size(I,2);


angstd=circshift(sawtooth([0:stimLen:simLen-stimLen]/1e3*2*pi*sacFreq,1)+1,[0 1]);
angstd=angstd*pi/24e3;
sMod=[.3 .7];


phiOffset=(rand(gridSz)-0.5)*pi/16;
phiOffset=phiOffset-mean(phiOffset(~logical(0+eye(gridSz))));
% phiOffset=phiOffset-diag(diag(phiOffset));

stimWeigth=ones(gridSz);
stimWeigth(ceil((gridSz+1)/2)+[-1: 1],ceil((gridSz+1)/2)+[-1: 1])=1;
% stimWeigth(ceil((gridSz+1)/2),ceil((gridSz+1)/2))=0;
stimWeigth(:,[1:gridSz]<=1+gridSz/2)=0;
stimWeigth(1,:)=1;


noiseWeigth=~stimWeigth*.0+.1;
noiseWeigth=noiseWeigth*5;

noiseStim=~(0*eye(gridSz))*pi/16;

for n=1:numel(phi)
  maxStr=(rand(gridSz)-.5)*.25+2;
    
  for x=1:gridSz
    selX=X==x;
    for y=1:gridSz
      phiDum=phi(n)+phiOffset(y,x)+randn*noiseStim(y,x);
      
      selY=Y==y;
      selXY=selX & selY & EI;
      
      for th=1:numOrient      
        I(selXY & theta==th, n)=stimWeigth(y,x)*.5*(cos(2*(th/numOrient*pi-phiDum))+1)*maxStr(y,x);
        I(selXY & theta==th, n)=I(selXY & theta==th, n)+noiseWeigth(y,x)*maxStr(y,x)*(rand(sum(selXY & theta==th),1)+.5);
      end
     
    end
  end
end

% I=.5*(I+fliplr(I));

I=I*1+(I>0)*2;

strMod=sMod(1)*(1-angstd/max(angstd(:)))+sMod(2);
I=bsxfun(@times,I,strMod);

snr=1;

noise=(I>0).*(I)/snr;

figure(9)
subaxis(2,2,1)
imagesc(tbins(1:end-1),1,I); colorbar
title('input')
subaxis(2,2,3)
imagesc(tbins(1:end-1),1,noise); colorbar
title('noise')
subaxis(3,2,2)
plot(tbins(1:end-1), mod(angle(exp(2i*phi))/2,pi))
title('input orientation')
subaxis(3,2,2,2)
plot(tbins(1:end-1), angstd)
title('standard deviation of orientation')
subaxis(3,2,2,3)
plot(tbins(1:end-1), strMod)
title('input strength')

inputDec=nan(gridSz,gridSz,size(I,2));
for n=1:size(I,2)
  inputDec(:,:,n)=decode(I(EI,n),theta(EI),X(EI),Y(EI),numOrient);
end




%%
% parameters RS
a=.02;
b=.2;
c=-65;
d=8;
parE=[a; b; c; d];

% parameters FS
a=.1;
b=.2;
c=-65;
d=2;
parI=[a; b; c; d];



numNeur=numel(EI);
V_init=repmat([c; b*c],1,numNeur);
V_init=V_init+bsxfun(@times,randn(1,numNeur),[1; .2]);

pars=nan(4,numNeur);
for n=1:4
  pars(n,EI)=parE(n)*(1+.01*randn(1,sum(EI)));
  pars(n,~EI)=parI(n)*(1+.01*randn(1,sum(~EI)));
end





tLim=minmax(tbins);
dt=.5;

input=[];
% input.I=[tbins(:) [I'; ]];
input.I=[tbins(:) [I'; I(:,1)']];
% input.noise = [tbins(:) [noise'; ]];
input.noise = [tbins(:) [noise'; noise(:,1)']];

input.EI=EI;





input.S=conMat;
input.STDP=0;
input.a=pars(1,:);
input.b=pars(2,:);
input.c=pars(3,:);
input.d=pars(4,:);
input.V_init=V_init;

input.tLim=tLim;
input.dt=dt;

input.EI=EI;

% input.fullOutput=1;


input.delay=0;

% input.saveInterval=1e4;

% input.output=['~/GIT/Dynamic_Coding/Temp_Bart/test_Orient/N40_8Col_8HypCol_2.mat'];

[output,spikes]=Izh_network_TAH(input);

% cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
% qsubfeval('Izh_network_TAH',input,'memreq',1024^3*6,'timreq',60^2*2);



%%
stimOnsets=[0:1e3/sacFreq:tLim(2)];
plotLim=[0 2]*1e3;

dt=output.input.dt;
EI=output.EI;

input=output.input;


ksig=2/dt;
tt=-5*ksig:5*ksig;
sphistkernel=exp(-tt.^2/(2*ksig^2));
sphistkernel=sphistkernel/sum(sphistkernel);

X=idx.X;
Y=idx.Y;

rateE=nan(gridSz,gridSz,numOrient,numel(output.t));
rateI=nan(gridSz,gridSz,numOrient,numel(output.t));
for x=1:gridSz
  selx=X==x;
  for y=1:gridSz
    selxy=selx & Y==y;
    for th=1:numOrient
      seltot=selxy & theta==th;
      rateE(y,x,th,:)=sum(output.spiks(EI & seltot,:));
      rateI(y,x,th,:)=sum(output.spiks(~EI & seltot,:));
    end
  end
end

figure(2)
set(2,'position',[100 100 1200 800]) 
clf
subaxis(5,1,1)
imagesc(output.t,1,conv2(reshape(squeeze(sum(rateE,3)),gridSz^2,[])',sphistkernel(:),'same').')
vline(stimOnsets)
ylabel('hypercolumn #')

subaxis(5,1,2)
imagesc(output.t,1,conv2(reshape(squeeze(sum(rateI,3)),gridSz^2,[])',sphistkernel(:),'same').')
xlabel('time (ms)')
xlim([output.input.tLim(1) max(output.t)])
xlim(plotLim)
vline(stimOnsets)
ylabel('hypercolumn #')

subaxis(5,1,3)
ns_plotspikes(spikes,gca, [], plotLim);
xlim([output.input.tLim(1) max(output.t)])
xlim(plotLim)
vline(stimOnsets)


% sphist=cat(1,full(sum(output.spiks(EI,:))),full(sum(output.spiks(~EI,:)))).';
[sp,f]= bg_fftwelch(reshape(squeeze(sum(rateE,3)),gridSz^2,[]).',1e3/dt,.5,1,1,1e3/dt);
subaxis(5,1,5)
imagesc(f(f>=0),1,mean(abs(sp(f>=0,:,:)).^2,3).')
xlim([0 100])
ylabel('hypercolumn #')
xlabel('freq (Hz)')

subaxis(5,1,4)
sphist=cat(1,full(sum(output.spiks(EI,:))),full(sum(output.spiks(~EI,:)))).';
plot(output.t,conv2(sphist,sphistkernel(:),'same'))


% subaxis(5,4,1,4,1,2)
% imagesc(reshape(mean(reshape(sum(output.spiks(EI,:),2)/((diff(output.I([1 end],1))+output.input.dt)/1e3),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
% colorbar
% axis image
% title('rate E')
% 
% subaxis(5,4,2,4,1,2)
% imagesc(reshape(mean(reshape(sum(output.spiks(~EI,:),2)/((diff(output.I([1 end],1))+output.input.dt)/1e3),cfg.numI,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
% colorbar
% axis image
% title('rate I')
% 
% dum=nan(cfg.gridSz,cfg.gridSz,numel(unique(idx.theta)));
% dumAng=dum;
% for n=1:numel(unique(idx.theta))
%   dum(:,:,n)=reshape(mean(reshape(sum(output.spiks(EI & idx.theta==n,:),2),[],cfg.gridSz^2)),cfg.gridSz,cfg.gridSz);
%   dumAng(:,:,n)=exp(1i*2*pi*(n/(numel(unique(idx.theta)))));
% end
% th_dec=sum(dum.*dumAng,3);
% 
% subaxis(5,4,3,4,1,2)
% imagesc(abs(th_dec))
% axis image
% colorbar
% 
% subaxis(5,4,4,4,1,2)
% imagesc(angle(th_dec)/2,[-pi pi]/2)
% axis image
% colorbar
% % colormap jet
% 
% figure(21)
% imagesc(angle(th_dec)/2,[-pi pi]/2)
% axis image
% colormap hsv
% %

figure(20)
clf
[ne,xe]=hist(sum(output.spiks(EI,:),2)/simLen*1e3,round(mean(sum(output.spiks(EI,:),2)/simLen*1e3))+[-5:.2:5]);
[ni,xi]=hist(sum(output.spiks(~EI,:),2)/simLen*1e3,round(mean(sum(output.spiks(~EI,:),2)/simLen*1e3))+[-10:.2:10]);
plot(xe,ne)
hold all
plot(xi,ni)
legend('E-neurons','I-neurons')
xlabel('firing rate (s^{-1})')
ylabel('N')

%% gamma phase coding

% determine gamma peak freq
[sp,f]=bg_fftwelch(reshape(rateE,[],size(rateE,4)).',1e3/dt,.5,1,1);
fsel=f>10 & f<80;
sp=sp(fsel,:,:);
f=f(fsel);
% go for common peak frequency
[~,midx]=max(mean(mean(abs(sp),3),2));
gamFreq=f(midx);

[tfr,~,t]=ft_specest_tfr(reshape(rateE,[],size(rateE,4)),output.t/1e3,'freqoi',gamFreq);
tfr=reshape(squeeze(tfr),size(rateE));

% gamma amplitude/power coding
amplDec=nan(gridSz,gridSz,size(tfr,4));
for x=1:gridSz
  for y=1:gridSz
    for n=1:size(tfr,4)      
      amplDec(x,y,n)=decode(squeeze(abs(tfr(x,y,:,n))),[1:8]',ones(8,1),ones(8,1),numOrient);
    end
  end
end

% phase coding
phasDec=nan(gridSz,gridSz,size(tfr,4));
for x=1:gridSz
  for y=1:gridSz
    for n=1:size(tfr,4)      
      decDum=squeeze(angle(tfr(x,y,:,n)./tfr(x,y,1,n)));
      decDum=decDum-min(decDum);
      phasDec(x,y,n)=decode(decDum,[1:8]',ones(8,1),ones(8,1),numOrient);
    end
  end
end
%% rate code
spiks=output.spiks(EI,:);
theta=idx.theta;
X=idx.X;
Y=idx.Y;
winLen=100/dt;
tSel=.5*winLen:10:tLim(2)/dt-.5*winLen;
winSel=[-winLen*0.5+1:winLen*0.5];
kernel=gaussKern(winLen/3,winLen,1)*1;
z=nan(gridSz,gridSz,numel(tSel));
for n=1:numel(tSel)
  selVec=tSel(n)+winSel;
  
  rate=sum(bsxfun(@times,spiks(:,selVec),kernel.'),2);
  
  z(:,:,n)=decode(rate,theta(EI),X(EI),Y(EI),numOrient);
  
  
end

tIntp=[tLim(1):dt:tLim(2)]';
phiInt=interp1(tbins(1:end-1)',exp(1i*2*phi(:)),tIntp);
decInt=interp1(tSel(:)*dt,squeeze(sum(sum(z,1),2)),tIntp);
decInt_centre=interp1(tSel(:)*dt,squeeze((z(3,3,:))),tIntp);
decCons=interp1(tSel(:)*dt,squeeze(abs(mean(mean(z,1),2)))./squeeze(mean(mean(abs(z),1),2)),tIntp);
decInp=interp1(tbins(1:end-1)',squeeze(sum(sum(inputDec,1),2)),tIntp);

figure(991)
clf
subaxis(2,1,1)
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(decInt./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
hold all
% plot(abs(nanmean((bg_reshape_overlap((0.5*angle(decInt_centre./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))),2)));
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(squeeze(mean(mean(amplDec(:,[1:gridSz]<=gridSz/2,:),1),2))./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(squeeze(mean(mean(phasDec(:,[1:gridSz]<=gridSz/2,:),1),2))./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
h=legend('rate code','gamma amplitude','gamma phase');
set(h,'location','best')
ylabel('\surd{\langle(\Delta\phi)^2\rangle}');
subaxis(2,1,2)
plot(tIntp,mod(0.5*angle(decInp),pi),'--')
hold all
plot(t*1e3,mod(0.5*angle(decInt),pi));
plot(t*1e3,mod(0.5*angle(squeeze(mean(mean(amplDec(:,[1:gridSz]<=gridSz/2,:),1),2))),pi));
plot(t*1e3,mod(0.5*angle(squeeze(mean(mean(phasDec(:,[1:gridSz]<=gridSz/2,:),1),2))),pi));
h=legend('stimulus','rate code','gamma amplitude','gamma phase');
set(h,'location','best')
ylabel('\phi');


figure(992)
plot([0:round(1000/sacFreq/dt)-1]*dt,(nanmean((bg_reshape_overlap((0.5*angle(squeeze(mean(mean(phasDec,1),2))./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))),2)));
hline(0,'-k')
ylabel('\langle\Delta\phi\rangle')

figure(993)
clf
hold all
plot(t*1e3,var(exp(1i*angle(rateDec)).'))
plot(t*1e3,var(exp(1i*angle(reshape(phasDec,16,[])))))
plot(t*1e3,var(exp(1i*angle(reshape(amplDec,16,[])))))
%%
sel=isnan(decInt) | isnan(decInt); 
[xc,lags]=(xcorr(exp(1i*decInt(~sel)),exp(1i*decInt(~sel)))); 
figure(13); 
plot(lags*dt,abs(xc))
[~,midx]=max(abs(xc));
vline(lags(midx)*dt);
title(num2str(lags(midx)*dt))

figure(12); 
clf

subaxis(4,2,1)
cla
hold all
try
% plot(tIntp,mod(0.5*angle(phiInt),pi),'--')
plot(tIntp,mod(0.5*angle(decInp),pi),'--')
% plot(tbins(1:end-1),fliplr(phi),'--')
end
plot(tIntp,mod(0.5*angle(decInt),pi),'k','linewidth',1)

vline(stimOnsets)

axis tight
legend('stim','decoder')
ylabel('\phi')
vline(stimOnsets)
ylim([0 pi])

subaxis(4,2,1,2)
plot(tIntp,(0.5*angle(decInt./decInp)))
hold all
plot(tIntp,(0.5*angle(decInt_centre./decInp)))
ylabel('\Delta \phi')
vline(stimOnsets)
hline(0,'-k')

subaxis(4,2,2,1,1,2)
plot(tSel*dt,squeeze(abs(mean(mean(z,1),2)))./squeeze(mean(mean(abs(z),1),2)))
hold on
% plot(tbins(1:end-1), angstd,'r')
title('consistency across hyper columns')
vline(stimOnsets)

drawnow

[sp,f,tt]=ft_specest_tfr(reshape(rateE,[],size(rateE,4)),[input.tLim(1):dt:input.tLim(2)]/1e3,'freqoi',[1:100]);
[sp2,f,tt]=ft_specest_tfr(sum(reshape(rateE,[],size(rateE,4))),[input.tLim(1):dt:input.tLim(2)]/1e3,'freqoi',[1:100]);

subaxis(4,1,3)
imagesc(tt*1e3,f,squeeze(abs(squeeze(sp2)))); axis xy
xlabel('time (ms')
vline(stimOnsets)

subaxis(4,1,4)
imagesc(tt*1e3,f,squeeze(mean(abs(squeeze(sp))))); axis xy
xlabel('time (ms')
vline(stimOnsets)

figure(14)
clf
subaxis(2,1,1)
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(decInt./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
hold all
% plot(abs(nanmean((bg_reshape_overlap((0.5*angle(decInt_centre./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))),2)));
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(squeeze(mean(mean(amplDec,1),2))./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
plot([0:round(1000/sacFreq/dt)-1]*dt,sqrt(nanmean((bg_reshape_overlap((0.5*angle(squeeze(mean(mean(phasDec,1),2))./decInp)),round(1000/sacFreq/dt),round(1000/sacFreq/dt))).^2,2)));
h=ylabel('\surd{\langle(\Delta\phi)^2\rangle}');
h=legend('rate code','gamma amplitude','gamma phase');
set(h,'location','best')

subaxis(2,1,2)
plot([0:round(1000/sacFreq/dt)-1]*dt,nanmean(bg_reshape_overlap(decCons,round(1000/sacFreq/dt),round(1000/sacFreq/dt)),2));
ylabel('consistency across hyper columns')
%%
figure(22)
for n=1:size(z,3)
  dum=z(:,:,n).';
  dum=sum(dum(:))/sum(abs(dum(:)));
%   for k=1:gridSz^2
%     subaxis(gridSz,gridSz,k)
%     compass(dum(k))
%   end
  h=compass(exp(1i));
  set(h,'visible','off')
  hold on
  
  compass(abs(dum)*exp(1i*angle(dum)/2))
  compass(-abs(dum)*exp(1i*angle(dum)/2))
  hold off
  title(num2str(tSel(n)*dt))
  drawnow
end