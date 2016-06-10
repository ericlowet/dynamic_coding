cfg=[];
cfg.numOrient=8;
cfg.gridSz=[1 1];
numE=20;
cfg.numE=numE;
% [nEE nEI nIE nII]

connFac=16;
sigFac=2;

cfg.numConIntra=[160 15 20 20]*connFac;
cfg.numConInter=[160 20 0 0]*connFac;
cfg.sig_dist=[2 2 .5 0];
cfg.sig_theta=[15 .15 .2 nan]*pi/sigFac;
% sigmas: intraE, intraI, interE, intrerI
cfg.strConIntra=[2e-3 2e-2 3e-2 3e-2]/connFac;
cfg.strConInter=[1e-3 2e-2 0 0]/connFac;
% cfg.strConInter=[1e-3 0 0 0];
cfg.noAutapse=0;

[conMat, idx, cfg, conMatSep] = genOrientColumns(cfg);


%% symmetric shifting
conMatComb=conMat;
shift=2;
for n=1:cfg.numOrient
  selE=true(cfg.numOrient*numE,1);
  selE((n-1)*numE+[1:numE])=false;
  conMatComb(selE,~selE)=circshift(conMat(selE,~selE),[shift,0]*numE);
%   conMatComb(~selE,selE)=circshift(conMat(~selE,selE),[0, shift]*numE);
end

%% creating hubs
conMatComb=conMat;
shiftSz=1;

centreNode=2;

for n=1:cfg.numOrient
  shift=sign(angle(exp(1i*centreNode*2*pi/8)*exp(-1i*n*2*pi/8)))*shiftSz;
  if round(angle(exp(1i*centreNode*2*pi/8)*exp(-1i*n*2*pi/8))*10)==round(-pi*10);
    shift=0;
  end
  selE=true(cfg.numOrient*numE,1);
  selE((n-1)*numE+[1:numE])=false;
  conMatComb(selE,~selE)=circshift(conMat(selE,~selE),[shift,0]*numE);
%   conMatComb(~selE,selE)=circshift(conMat(~selE,selE),[0, shift]*numE);
end



figure(8)
subaxis(2,2,1)
imagesc(conMat>0)
title('complete S')
subaxis(2,2,1,2)
imagesc(conMatComb>0)
title('complete S')


subaxis(2,2,2)
imagesc(log(conMat))
% imagesc(conMat,[0 quantile(conMat(conMat(:)>0),.9)])
title('complete S')
axis image
h=colorbar;
set(get(h,'ylabel'),'string','Log(S)')

subaxis(2,2,2,2)
imagesc(log(conMatComb))
% imagesc(conMat,[0 quantile(conMat(conMat(:)>0),.9)])
title('complete S')
axis image
h=colorbar;
set(get(h,'ylabel'),'string','Log(S)')

%%
savcfg=1;
clear basedir
if savcfg
  basedir=['~/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/29/'];
%   basedir=['~/GIT/Dynamic_Coding/Temp_Bart/test_Orient/N' num2str(numE) num2str([sigFac connFac],'_sigFVar_connF%d')];
%   basedir=['~/GIT/Dynamic_Coding/Temp_Bart/test_Orient/N' num2str(numE) num2str([sigFac connFac],'_sigF%d_connF%d')];
  if ~exist(basedir,'file')
    mkdir(basedir)
  end
  if ~exist(fullfile(basedir,'stim'))
    mkdir(fullfile(basedir,'stim'))
  end
  cfg.S=conMat;
  save(fullfile(basedir,'stim','aux.mat'),'idx','cfg')
end
%%
X=idx.X;
Y=idx.Y;
EI=idx.EI;
theta=idx.theta;
numNeur=numel(X);
numOrient=cfg.numOrient;
gridSz=cfg.gridSz;


cfgms=[];
cfgms.sampling_frequency=1e3;
cfgms.ms_interval=400;
cfgms.simulation_time=cfgms.ms_interval*8*5;
cfgms.ms_interval_var=0;
cfgms.modulation_strength=.75;
cfgms.modulation_var=1e-2;
cfgms.positive_mod=4/3;
cfgms.negative_mod=3/3;
cfgms.base_rate=1;
[ ms_signal,ms_times]  = create_ms_sig(cfgms);
ms_signal=max(circshift(ms_signal,18),0);

sacFreq=1000/cfgms.ms_interval;
isi=cfgms.ms_interval;
% stimLen=10;
numStim=floor(cfgms.simulation_time/isi);
simLen=isi*numStim;


% tuning functions
figure(7)
thh=linspace(-pi/2,pi/2,1001);
th2=([1:numOrient]/numOrient*2*pi-pi)/2;

TF=(cos(thh*2)+1)/2;
tf_sig=pi/8;
TF=exp(-.5*thh.^2/tf_sig^2);
TF2=exp(-.5*th2.^2/tf_sig^2);
plot(thh,TF, th2, TF2)
title('Tuning function')
ylabel('input strength')
xlabel('stimulus orientation')
axis tight

stimOrient=[1:5].';
stimOrient=stimOrient*2*pi/numel(stimOrient);
phiIdxdum=repmat([1:numel(stimOrient)],numStim/numel(stimOrient),1);
phiIdx=nan(numel(phiIdxdum),prod(gridSz));
phiIdx(:,1)=phiIdxdum(:);
if size(phiIdx,2)==2
phiIdxdum=phiIdxdum';
phiIdx(:,2)=phiIdxdum(:);
end
% shift second stimulus 90 degrees
% phiIdx(:,2)=mod(phiIdx(:,2)+3,8)+1;

%randomly permute the stimuli
phiIdx=phiIdx(randperm(numStim),:);

phi=stimOrient(phiIdx);


tdum=[1:isi:simLen; isi:isi:simLen];

tbins=min(tdum(:)):1e3/cfgms.sampling_frequency:max(tdum(:));


I=zeros(numNeur,numel(tbins));
I(EI,:)=1;




noiseWeigth=ones(gridSz)*.2;
baseLine=ones(gridSz)*5;
stimWeigth=ones(gridSz)*1;
% stimWeigth=[1.3; 0.7];

for n=1:size(tdum,2)
  for x=1:gridSz(1)
    selX=X==x;
    for y=1:gridSz(2)
      phiDum=phi(n,x,y);%+randn*noiseStim(x,y);
      
      selY=Y==y;
      selXY=selX & selY & EI;
      
      for th=1:numOrient      
        dum=stimWeigth(x,y)*exp(-.5*(.5*angle(exp(1i*(2*th/numOrient*pi-phiDum)))).^2/tf_sig^2);
        [~,midx1]=min(abs(tbins-tdum(1,n)));
        if n<size(tdum,2)
          [~,midx2]=min(abs(tbins-tdum(2,n+1)));
        else
          midx2=max(tbins);
        end
        selVec=midx1:midx2;
        dum=dum+baseLine(x,y)+noiseWeigth(x,y)*(rand(sum(selXY & theta==th),numel(selVec))+.5);
        I(selXY & theta==th, selVec)=dum;
      end
     
    end
  end
%   I(~EI,n)=(rand(sum(~EI),1)-.5)*.25+1;
end


% % evoked response shape
% tt=[0:stimLen:isi];
% 
% ka=30;
% t0a=0;
% alph_freq=(3-2./(1+exp(-ka*(tt-t0a))));
% 
% tau=[25 50];
% A=-(tau(2)/tau(1))^(-tau(2)/(tau(2)-tau(1)))+(tau(2)/tau(1))^(-tau(1)/(tau(2)-tau(1)));
% alph_env=1/A*(exp(-tt/tau(2))-exp(-tt/tau(1)));
% alph_env(tt<0)=0;
% 
% tau=[25 100];
% A=-(tau(2)/tau(1))^(-tau(2)/(tau(2)-tau(1)))+(tau(2)/tau(1))^(-tau(1)/(tau(2)-tau(1)));
% inp_env=1/A*(exp(-tt/tau(2))-exp(-tt/tau(1)));
% inp_env(tt<0)=0;
% 
% kg=10;
% t0g=.2;
% 
% strModDum=nan(numel(tt),numStim);
% for m=1:numStim
% modulation=sin(tt/1e3*2*pi*(10+randn*1.5).*alph_freq).*alph_env;
% strModDum(:,m)=0.5+modulation+inp_env;
% end
% strMod=strModDum(1:end-1,:);

strMod=ms_signal;
I=bsxfun(@times,I,strMod(:)');

snr=.5;

noise=(I>0).*(I)/snr;

figure(9)
subaxis(2,2,1)
imagesc(tbins(1:end-1),1,I); colorbar
title('input')
xlim([0 min(8*cfgms.ms_interval, max(tbins))])
subaxis(2,2,3)
imagesc(tbins(1:end-1),1,noise); colorbar
title('noise')
xlim([0 min(8*cfgms.ms_interval, max(tbins))])
subaxis(2,2,2)
plot(tdum(1,:), phi/pi*180,'.-')
title('input orientation')
axis tight
ylim([0 360])
set(gca,'ytick',[0:45:360])
% subaxis(3,2,2,2)
% plot(tbins(1:end), angstd)
% title('standard deviation of orientation')
subaxis(2,2,2,2)
plot(tbins(1:end), strMod(:))
title('input strength')
axis tight
xlim([0 cfgms.ms_interval*8])

inputDec=nan(gridSz(1),gridSz(2),size(I,2));
for n=1:size(I,2)
  inputDec(:,:,n)=decode(I(EI,n),theta(EI),X(EI),Y(EI),numOrient);
end
dt=.5;
tLim=minmax(tbins);
tIntp=[tLim(1):dt:tLim(2)]';

decInp=interp1(tbins(1:end)',squeeze(inputDec).',tIntp);
%%
savStim=1;
if savStim
  if ~exist(fullfile(basedir,'stim'))
    mkdir(fullfile(basedir,'stim'))
  end
  save(fullfile(basedir,'stim','decInp.mat'),'decInp','phiIdx','cfgms')
end

return
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


input=[];
% input.I=[tbins(:) [I'; ]];
input.I=[tbins(:) [I'; ]];
% input.noise = [tbins(:) [noise'; ]];
input.noise = [tbins(:) [noise'; ]];

input.EI=EI;






input.STDP=0;
input.a=pars(1,:);
input.b=pars(2,:);
input.c=pars(3,:);
input.d=pars(4,:);
input.V_init=V_init;

input.tLim=tLim;
input.dt=dt;

input.EI=EI;
input.S=conMat;
 
input.fullOutput=0;
input.LFPout=1;
input.LFPkernel=[];
for n=1:cfg.numOrient
  input.LFPkernel=[input.LFPkernel (EI & idx.theta==n)];
end

input.S=conMat;
input.delay=0;

input.saveInterval=1e4;
%%
shifts=[-cfg.numE:cfg.numE/2:cfg.numE];
shifts=-2:2;
shifts=1:5;
numConnPat=numel(shifts);
nN=numE*cfg.numOrient;
for nn=1:numConnPat
  input.output=fullfile(basedir, ['N' num2str(cfg.numE) '_' num2str(cfg.numOrient) 'Col_P' num2str(nn,'%02d') '.mat']);
  if ~exist(input.output,'file')
    
    centreNode=shifts(nn);
    
  
    
    conMatComb=conMat;
    shift=shifts(nn);
    for n=1:cfg.numOrient
      
      % when doing asymmetric shifting
      shift=sign(angle(exp(1i*centreNode*2*pi/8)*exp(-1i*n*2*pi/8)))*shiftSz;
      if round(angle(exp(1i*centreNode*2*pi/8)*exp(-1i*n*2*pi/8))*10)==round(-pi*10);
        shift=0;
      end
      %
      
      selE=true(cfg.numOrient*numE,1);
      selE((n-1)*numE+[1:numE])=false;
      conMatComb(selE,~selE)=circshift(conMat(selE,~selE),[shift,0]*numE);
    end
    
    
    
    input.S=conMatComb;
    %     input.S=conMat;
    %     input.S(EI,EI)=conMatComb;
    
    
%       [output,spikes]=Izh_network_TAH(input);
    
    input.verbose=0;
    cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
    qsubfeval('torqueWrap','Izh_network_TAH',input,'memreq',1024^3*2,'timreq',60*60);
    
  end
end
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

rateE=nan(gridSz(1),gridSz(2),numOrient,numel(output.t));
rateI=nan(gridSz(1),gridSz(2),numOrient,numel(output.t));
for x=1:gridSz(1)
  selx=X==x;
  for y=1:gridSz(2)
    selxy=selx & Y==y;
    for th=1:numOrient
      seltot=selxy & theta==th;
      rateE(x,y,th,:)=sum(output.spiks(EI & seltot,:));
      rateI(x,y,th,:)=sum(output.spiks(~EI & seltot,:));
    end
  end
end

figure(2)
set(2,'position',[100 100 1200 800]) 
clf
subaxis(5,1,1)
imagesc(output.t,1,conv2(reshape(squeeze(sum(rateE,3)),prod(gridSz),[])',sphistkernel(:),'same').')
vline(stimOnsets)
ylabel('hypercolumn #')

subaxis(5,1,2)
imagesc(output.t,1,conv2(reshape(squeeze(sum(rateI,3)),prod(gridSz),[])',sphistkernel(:),'same').')
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
[spfft,f]= bg_fftwelch(reshape(squeeze(sum(rateE,3)),prod(gridSz),[]).',1e3/dt,.5,1,1,1e3/dt);
subaxis(5,1,5)
plot(f(f>=0),mean(abs(spfft(f>=0,:,:)).^2,3).')
xlim([0 100])
% ylabel('hypercolumn #')
ylabel('power')
xlabel('freq (Hz)')

subaxis(5,1,4)
sphist=cat(1,full(sum(output.spiks(EI,:))),full(sum(output.spiks(~EI,:)))).';
plot(output.t,conv2(sphist,sphistkernel(:),'same'))

figure(20)
clf
simLen1=diff(minmax(output.t));
[ne,xe]=hist(sum(output.spiks(EI,:),2)/simLen1*1e3,round(mean(sum(output.spiks(EI,:),2)/simLen1*1e3))+[-5:.2:5]);
[ni,xi]=hist(sum(output.spiks(~EI,:),2)/simLen1*1e3,round(mean(sum(output.spiks(~EI,:),2)/simLen1*1e3))+[-10:.2:10]);
plot(xe,ne)
hold all
plot(xi,ni)
legend('E-neurons','I-neurons')
xlabel('firing rate (s^{-1})')
ylabel('N')


[sp,f,ttt]=ft_specest_tfr(squeeze(sum(rateE,3)),output.t/1e3,'freqoi',1:120);
figure(21)
subaxis(2,1,1)
imagesc([1:isi],f,abs(nanmean(reshape(squeeze((sp(1,:,1:(isi/dt)*floor(numel(ttt)/(isi/dt))))),numel(f),isi/dt,[]),3)))
axis xy
subaxis(2,1,2)
imagesc([1:isi],f,nanmean(reshape(squeeze(abs(sp(1,:,1:(isi/dt)*floor(numel(ttt)/(isi/dt))))),numel(f),isi/dt,[]),3))
% imagesc(ttt,f,squeeze(abs(sp(2,:,:))))
axis xy

