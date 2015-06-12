addpath ~/GIT/Dynamic_Coding/Temp_Bart/paruly





numNeurLin=50;

dR=50e-6;
% E->E, E->I, I->E, I->I
sig=[3 10 3 3]*dR;
p=[.2 .5 .8 .5];
strength=[.1 1 1 1];
[ connMat, connMat_sep, R_e, R_i ] = genConnPING(numNeurLin, dR, sig, p, strength);
numE=size(R_e,1);
numI=size(R_i,1);
S=connMat;
% S(1:numE,1:numE)=double(connMat_sep{1}>0)*.01;

%% generate stimuli
saveStim=0;
genStim=0;
if genStim
  numStim=5;
  corrLim=0;
  I=randn(numNeurLin,numNeurLin,numStim);
  sigI=4*dR;
  Rdum=min(R_e):dR:max(R_e);
  Rdum=Rdum-mean(Rdum);
  filtKern=exp(-Rdum.^2/(2*sigI^2));
  filtKern=bsxfun(@times,filtKern,filtKern');
  filtKern=filtKern/sum(filtKern(:));
  for n=1:numStim
    Idum=I(:,:,n);
    Idum=ifft2(fft2(Idum).*fft2(filtKern));
    Idum=Idum-min(Idum(:));
    Idum=Idum/mean(Idum(:))*4;
    I(:,:,n)=Idum;
  end
  
  Idum=permute(I,[3 1 2]);
  R=corr(Idum(:,:)');
  R(logical(eye(numStim)))=nan;
  mR=nanmax(R);
  
  [mR,mIdx]=max(mR);
  genFlag= mR > corrLim;
  while genFlag
    
    Idum=rand(numNeurLin);
    Idum=ifft2(fft2(Idum).*fft2(filtKern));
    Idum=Idum-min(Idum(:));
    Idum=Idum/mean(Idum(:))*4;
    I(:,:,mIdx)=Idum;
    
    Idum=permute(I,[3 1 2]);
    R=corr(Idum(:,:)');
    R(logical(eye(numStim)))=nan;
    mR=nanmax(R);
    [mR,mIdx]=max(mR);
    genFlag= mR > corrLim;
  end
  
  Idum=permute(I,[3 1 2]);
  R=corr(Idum(:,:)');  
  
  stims=I;
  %% generate blocks
  blockLen=1e3; %1 second blocks
  dt=.5;
  tLim=[0 50e3]; % 50 seconds, --> 10 seconds per block
  
  numBlocks=round(tLim(2)/blockLen);
  
  numTrains=ceil(numBlocks/numStim);
  
  trains=repmat([1:5]',numTrains,1);
  trains=trains(1:numBlocks);
  
  trains=[trains'; trains'];
  trains=trains(:); %double up: once for onset, once for ofset. Otherwise the simulation would interpolate between stimuli
  
  tim=[dt:blockLen:tLim(2); blockLen:blockLen:tLim(2)];
  
  Idum=permute(stims,[3 1 2]);
  
  I=zeros(numTrains*2,numE+numI);
  
  for n=1:numel(trains)
    I(n,1:numE)=Idum(trains(n),:);
  end
  
  % add  timing column
  I=[tim(:) I];
  
  % SNR = 10
  noise=I.*.01;
  noise(:,1)=tim(:);
else
  numStim=5;
  fname=['N' num2str(numNeurLin) '_S' num2str(numStim)];
  load(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/stimuli',fname));
  Idum=permute(stims,[3 1 2]);
  R=corr(Idum(:,:)');
end


figure(1)
clf
set(1,'position',[100 100 1200 800])
colormap(paruly(128))
for n=1:numStim
subaxis(2,3,n)
imagesc(stims(:,:,n))
axis image
h=colorbar;
set(h,'location','southoutside')
title(['stimulus ' char(64+n)])
end
subaxis(2,3,n+1)

imageNan(1,1,R,[-1 1])
set(gca,'xtick',1:numStim,'xticklabel',char(64+[1:numStim]'))
set(gca,'ytick',1:numStim,'yticklabel',char(64+[1:numStim]'))
axis image
h=colorbar;
set(h,'location','southoutside')
title('Correlation in inputs')


if saveStim
  fname=['N' num2str(numNeurLin) '_S' num2str(numStim)];
  save(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/stimuli',fname),'stims','I','trains')
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/stimuli',fname),'-png','-eps','-transparent',1)  
end











%% neuron and simulation parameters
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



numNeur=size(R_e,1)+size(R_i,1);
V_init=repmat([c; b*c],1,numNeur);
V_init(1,:)=V_init(1,:)+randn(1,numNeur)*10;


numNeur=numE+numI;
EI=false(1,numNeur);
EI(1:numE)=true;

pars=[repmat(parE,1,numE) repmat(parI,1,numI)];
pars(3,1:numE)=pars(3,1:numE)+randn(1,numE)*3;
pars(4,1:numE)=pars(4,1:numE)+randn(1,numE)*1;
pars(1,numE+1:end)=pars(1,numE+1:end)+rand(1,numI)*.03;
pars(2,numE+1:end)=pars(2,numE+1:end)+rand(1,numI)*.03;

% plot synapstic connections
figure(2)
clf
colormap(paruly(128))
set(2,'position',[100 100 1000 500])
subaxis(1,2,1)
imagesc(S>0)
axis image
h=colorbar;
set(h,'location','southoutside')
title('adjancency matrix')
subaxis(1,2,2)
imagesc(S,[0 quantile(S(S>0),.9)])
axis image
h=colorbar;
set(h,'location','southoutside')
title('S (synaptic conductances)')


%%
cfg=[];
continued=0;
if continued
  load(fullfile('~/GIT/Dynamic_Coding/Temp_Bart/','N30_S5_snr0.1_50s_block_1s.mat'))
  cfg.S=output.S;
else
  cfg.S=S;
end


cfg.tLim=tLim;
cfg.dt=dt;
cfg.a=pars(1,:);
cfg.b=pars(2,:);
cfg.c=pars(3,:);
cfg.d=pars(4,:);

cfg.EI=EI;
cfg.STDP=1;
cfg.output=['/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/N' num2str(numNeurLin) '_S' num2str(numStim) '_snr0.1_' num2str(tLim(2)/1e3) 's_block_' num2str(blockLen/1e3) 's.mat'];
cfg.saveInterval=2.5e3;
cfg.verbose=0;
cfg.fullOutput=0;

cfg.A_STDP=[1 1]*1;


% tic
% [V,t,output,spikes]=Izh_network_STDP_izh2004(V_init,tLim,dt,I,noise,cfg);
% toc

% tic
% [output,spikes]=Izh_network(V_init,I,noise,cfg);
% toc

cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
qsubfeval('Izh_network',V_init,I,noise,cfg,'memreq',1024^3*12,'timreq',60^2*4);


%%
EI=output.EI;
spiks=output.spiks;
I=output.I;

try t=output.t; end
[spec,f,tt]=ft_specest_tfr(full(sum(spiks(EI,:))),t/1e3,'freqoi',[0:100],'verbose',0,'width',21);
  
tLimPlot2=[-.5e3 0]+max([spikes.timestamp{:}]);
tLimPlot1=[0 500]+2250;
dt=diff(t(1:2));

numE=sum(output.EI);
numI=numel(output.EI)-numE;

figure(10)
clf
set(10,'position',[100 100 1200 800])
ns_plotspikes(spikes,10,421);
xlim(tLimPlot1)
ns_plotspikes(spikes,10,422);
xlim(tLimPlot2)

kernel=exp(-[(-20/dt):(20/dt)].^2/(2*(1/dt)^2)); 
kernel=kernel/sum(kernel);


% subaxis(5,1,1,1,1,3)
% imagesc(t,1,-spiks)
% colormap gray
% xlim(tLimPlot)
% subaxis(5,1,4)
subaxis(4,2,1,2)
plot(t,conv(full(sum(spiks(~EI,:))),kernel,'same'),'color',[0 .5 0])
hold on
plot(t,conv(full(sum(spiks(EI,:))),kernel,'same'))
axis tight
xlim(tLimPlot1)
ylabel('# of neurons firing')

subaxis(4,2,2,2)
plot(t,conv(full(sum(spiks(~EI,:))),kernel,'same'),'color',[0 .5 0])
hold on
plot(t,conv(full(sum(spiks(EI,:))),kernel,'same'))
axis tight
xlim(tLimPlot2)


subaxis(4,1,3)
cla
% plot(f(f>=0),abs(spec(f>=0,:)))
% xlim([0 100])
imageNan(tt*1e3,f,squeeze(abs(spec)))
xlim([0 max([spikes.timestamp{:}])])
axis xy
% hold on
% boxCarLen=.3;
% h=plot(tt,conv(full(sum(spiks(EI,:))),ones(1,1e3*boxCarLen/dt),'same')/numE/boxCarLen,'color',[1 1 1]*.6,'linewidth',2);
% legend(h,'E firing rate')

% 
% subaxis(4,2,1,4)
% cla
% imagesc(reshape(output.I(1,[1:numE]+1),[sqrt(numE),sqrt(numE)]))


subaxis(4,2,1,4)
x=0:.5:80;
tSel=100:find(t==max([spikes.timestamp{:}]));
[Ne]=hist(sum(spiks(EI,tSel),2)/(numel(t(tSel))*dt*1e-3),x);
[Ni]=hist(sum(spiks(~EI,tSel),2)/(numel(t(tSel))*dt*1e-3),x);
% hold on
% bar(x,Ne/numE)
% bar(x,Ni/numI,'facecolor',[0 .5 0])
% xlim([3 max(x)])
% legend('E','I')
plot(x,[Ne(:)/numE,Ni(:)/numI])
xlabel('Firing rate (Hz)')
xlim([3 max(x)])
legend('E','I')
xlabel('Firing rate')
title('histogram of F-rates')

subaxis(4,2,2,4)
plot(output.I(1,2:numE+1),sum(output.spiks(1:numE,:),2)*1e3/max([spikes.timestamp{:}]),'.')
xlabel('input (nA)')
ylabel('Firing rate (Hz)')

fsz=14;
S_out=output.S;
S_orig=output.S_orig;
figure(11)
clf
set(11,'position',[150 150 1400 600])
subaxis(1,3,1,'mr',.02,'ml',.02)
imagesc(S_orig,[0, quantile(S_orig(S_orig>0),.9)])
axis image
h=colorbar;
set(h,'location','southoutside')
title('original S','fontsize',fsz)
subaxis(1,3,2)
imagesc(S_out,[0, quantile(S_orig(S_orig>0),.9)])
axis image
h=colorbar;
set(h,'location','southoutside')
title('new S','fontsize',fsz)
subaxis(1,3,3)
dS=S_out-S_orig;
% dS=dS(1:numE,1:numE);
imagesc(dS,max(abs(quantile(dS(abs(dS)>0),[.1 .9])))*[-1 1])
axis image
h=colorbar;
set(h,'location','southoutside')
title('\DeltaS','fontsize',fsz)
bg_redvbluecmap
% colormap([linspace(1,0,256)',ones(256,1),linspace(0,1,256)'])

savefig=0;
if savefig
  fname='scholarp_multip_norm_10000A'
%   fname='Izh_STDP_hard_bound';
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs',[fname '_output']),'-png','-pdf','-transparent',10);
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs',[fname '_Synap']),'-png','-pdf','-transparent',11);
end

%% generating decoding templates
numNeurLin=30;
numE=numNeurLin^2;
numStim=5;
blockLen=1e3; %1 second blocks
dt=.5;
tLim=[0 50e3]; % 50 seconds, --> 10 seconds per block
load(fullfile('~/GIT/Dynamic_Coding/Temp_Bart/',['N' num2str(numNeurLin) '_S' num2str(numStim) '_snr0.1_' num2str(tLim(2)/1e3) 's_block_' num2str(blockLen/1e3) 's_part2.mat']))
numStim=5;
fname=['N' num2str(numNeurLin) '_S' num2str(numStim)];
load(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/stimuli',fname));

numBlocks=round(tLim(2)/blockLen);

winLen=20/dt; %20 ms;
winLenT=numel([ceil(-winLen/2):floor(winLen/2)]);
dat=cell(numStim,1);
ydum=dat;

x=[1:numNeurLin]-(numNeurLin+1)/2;
y=x';
tt=permute([1:winLenT]-(winLenT+1)/2,[1 3 2]);
sig=[1 1 2];
kernel=exp(-x.^2/(2*sig(1)^2));
kernel=bsxfun(@times,kernel,exp(-y.^2/(2*sig(2)^2)));
kernel=bsxfun(@times,kernel,exp(-tt.^2/(2*sig(3)^2)));
kernel=kernel./sum(kernel(:));

uSamp=[1, 1];
for n=1:numStim
  tims=output.I(trains==n,1)/dt;
  tims=reshape(tims,2,[]);
  datExtr=cell(size(tims,2),1);
  for k=1:size(tims,2)
    datdum=full(output.spiks(1:numE,tims(1,k):tims(2,k)));
    [~,pidx]=findpeaks(sum(datdum),'minpeakdistance',15/dt);
%     datExtr{k}=zeros(ceil(1/4*numNeurLin)^2,ceil(winLenT/6),numel(pidx));    
    datExtr{k}=zeros(ceil(1/uSamp(1)*numNeurLin)^2,ceil(winLenT/uSamp(2)),numel(pidx)); 
    for p=1:numel(pidx)
      tSel=pidx(p)+[ceil(-winLen/2):floor(winLen/2)];
      tSelS=tSel>0 & tSel<size(datdum,2);
      datGam=zeros(numE,winLenT);
      datGam(:,tSelS)=datdum(:,tSel(tSelS));
      datGam=reshape(datGam,[numNeurLin numNeurLin winLenT]);
      
      %decrease feature space (smoothing and under sampling)
      datGamSmth=ifftn(fftn(datGam).*fftn(ifftshift(kernel)));
      datGamSmth=datGamSmth(1:uSamp(1):end,1:uSamp(1):end,1:uSamp(2):end);
      datExtr{k}(:,:,p)=reshape(datGamSmth,size(datExtr{k},1),[]);
      
    end
  end
  dat{n}=permute(cat(3,datExtr{:}),[3 1 2]);
  dat{n}=dat{n}(:,:);
	ydum{n}=ones(size(dat{n},1),1)*n;
end

template=nan(numE*winLenT,numStim);
for n=1:numStim
%   template(:,:,n)=reshape(mean(dat{n}),[numE,winLenT]);
  template(:,n)=mean(dat{n});
  template(:,n)=zscore(template(:,n));
end

%gram-schmidt orthonormalization
templateGS=template;
for n=1:numStim  
  templateGS(:,n)=templateGS(:,n)./norm(templateGS(:,n));  
  for k=n+1:numStim
    templateGS(:,k)=templateGS(:,k)-project(templateGS(:,k),templateGS(:,n));
  end  
end

%% showing templates

figure(4)
set(gcf,'position',[100 100 1200 900])
for n=1:numStim
  subaxis(4,5,n,1)
  imagesc(stims(:,:,n))
  if n==1
    ylabel('stimulus')
  end
  title(char(64+n))
  subaxis(4,5,n,2)
  imagesc(mean(reshape(template(:,n),[numNeurLin,numNeurLin,winLenT]),3))
  if n==1
    ylabel('spatial decoding template')
  end
  [mpk,midx]=max(reshape(template(:,n),[numNeurLin*numNeurLin,winLenT]),[],2);
  midx(mpk<0)=nan;
  subaxis(4,5,n,3)
  imageNan(reshape(midx,[numNeurLin,numNeurLin]),[minmax(midx(:))])
  if n==1
    ylabel('gamma spike timing')
  end
  subaxis(4,5,n,4)
  imagesc(reshape(template(:,n),[numNeurLin*numNeurLin,winLenT]))
  if n==1
    ylabel('spatio-temporal decoding template')
  end
end

saveFig=0;
if saveFig
  figdir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs/decoding';
  fname='templates';
  export_fig(fullfile(figdir,fname),'-eps','-png','-transparent',4)
end

%%
stim=ones(size(spiks,2),1);

for n=2:numStim
  for k=1:numBlocks/numStim
    tSel=[1:blockLen/dt]+(n-1)*blockLen/dt+(k-1)*numStim*blockLen/dt;
    stim(tSel)=n;
  end
end


%% decoding
spiks=output.spiks;
bulkLen=1e3;
numIts=ceil(size(spiks,2)/bulkLen);

projVec=nan(size(spiks,2),numStim,1);

for k=1:numIts
  disp([num2str(k) '/' num2str(numIts)])
  try
    a=bg_reshape_overlap(spiks(1:numE,[1:bulkLen+winLenT-1]+(k-1)*bulkLen),1,winLenT,2);
  catch
    a=bg_reshape_overlap(spiks(1:numE,((k-1)*bulkLen+1):end),1,winLenT,2);
  end
  dum=reshape(a,numE*winLenT,size(a,3));
  
  projVec([1:size(dum,2)]+(k-1)*bulkLen,:,1)=dum.'*template;
%   projVec([1:size(dum,2)]+(k-1)*bulkLen,:,2)=dum.'*templateGS;  
end
  
fname=['N' num2str(numNeurLin) '_S' num2str(numStim)];
load(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/stimuli',fname));
[~,decVec]=max(projVec,[],2);
decVec=squeeze(decVec);
accur=mean(bsxfun(@eq,decVec,stim));
%%


xlims=[0 1e3];
ctims=[1.05e4 3.75e4];
figure(3)  
clf
set(3,'position',[100 100 1200 800])
ns_plotspikes(spikes,3,[321])
xlim(xlims+ctims(1))
ns_plotspikes(spikes,3,[322])
xlim(xlims+ctims(2))
subaxis(3,2,1,2)
plot(output.t,projVec(:,:,1))
xlim(xlims+ctims(1))
subaxis(3,2,2,2)
plot(output.t,projVec(:,:,1))
xlim(xlims+ctims(2))
subaxis(3,2,1,3)
cla
plot(output.t,[decVec(:,1)])% decVec(:,2)])
h=legend(num2str(accur'));
set(get(h,'title'),'string','accuracy')
set(h,'location','northwest')
hold on
plot(output.t,stim,'--k','linewidth',2)
xlim(xlims+ctims(1))
xlabel('time (ms)')
subaxis(3,2,2,3)
plot(output.t,[decVec(:,1)]);% decVec(:,2)])
hold on
plot(output.t,stim,'--k','linewidth',2)
xlim(xlims+ctims(2))
xlabel('time (ms)')

saveFig=0;
if saveFig
  figdir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs/decoding';
  fname='performance_part1';
  export_fig(fullfile(figdir,fname),'-eps','-png','-transparent',3)
end

