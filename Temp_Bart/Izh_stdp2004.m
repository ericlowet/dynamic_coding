addpath ~/GIT/Dynamic_Coding/Temp_Bart/paruly
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

tLim=[0 10e3];
dt=.5;


numNeur=20;
dR=50e-6;
% E->E, E->I, I->E, I->I
sig=[2 3 10 3]*dR;
sig=[3 20 5 10]*dR;
p=[.1 .5 .8 .5];
strength=[.1 1 1 1];
[ connMat, connMat_sep, R_e, R_i ] = genConnPING(numNeur, dR, sig, p, strength);
numE=size(R_e,1);
numI=size(R_i,1);
S=connMat;
% S(1:numE,1:numE)=double(connMat_sep{1}>0)*.01;

% S=S_out;

I=randn(numNeur,numNeur);
sigI=2*dR;
Rdum=min(R_e):dR:max(R_e);
Rdum=Rdum-mean(Rdum);
filtKern=exp(-Rdum.^2/(2*sigI^2));
filtKern=bsxfun(@times,filtKern,filtKern');
filtKern=filtKern/sum(filtKern(:));
I=ifft2(fft2(I).*fft2(filtKern));
I=I-min(I(:));
I=I/mean(I(:))*3;

I_orig=I;

figure(1)
clf
set(1,'position',[100 100 1200 500])
colormap(paruly(128))
subaxis(1,3,1)
imagesc(I_orig)
axis image
h=colorbar;
set(h,'location','southoutside')
title('Input')
subaxis(1,3,2)
imagesc(S>0)
axis image
h=colorbar;
set(h,'location','southoutside')
title('adjancency matrix')
subaxis(1,3,3)
imagesc(S)
axis image
h=colorbar;
set(h,'location','southoutside')
title('S (synaptic conductances)')

savefig=0;
if savefig
  fname='inputs'
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs',[fname]),'-png','-pdf','-transparent',1);
end



I=[I(:); zeros(size(R_i,1),1)]';

I=[tLim(1) I; tLim(2) I];

noise=logical(I)*0;
noise(:,1)=tLim;


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


%%



cfg=[];
cfg.a=pars(1,:);
cfg.b=pars(2,:);
cfg.c=pars(3,:);
cfg.d=pars(4,:);
cfg.S=S;
cfg.EI=EI;
cfg.STDP=1;
% cfg.output='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test2_STDP_2004.mat';
cfg.saveInterval=inf;
cfg.verbose=1;

cfg.A_STDP=[1 1]*0.01;

% 
% 
tic
[V,t,output,spikes]=Izh_network_STDP_izh2004(V_init,tLim,dt,I,noise,cfg);
toc

tic
[V,t,output,spikes]=Izh_network(V_init,tLim,dt,I,noise,cfg);
toc

% cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
% qsubfeval('Izh_network',V_init,tLim,dt,I,noise,cfg,'memreq',1024^3*6,'timreq',60^2*2);


%%
EI=output.EI;
spiks=output.spiks;
I=output.I;
[spec,f,tt]=ft_specest_tfr(full(sum(spiks(EI,:))),t/1e3,'freqoi',[0:100],'verbose',0,'width',21);
  
tLimPlot=[-1e3 0]+tLim(2);
figure(10)
clf
set(10,'position',[100 100 1200 800])
ns_plotspikes(spikes,10,421);
xlim([0 1e3])
ns_plotspikes(spikes,10,422);
xlim(tLimPlot)

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
xlim([0 1e3])
ylabel('# of neurons firing')

subaxis(4,2,2,2)
plot(t,conv(full(sum(spiks(~EI,:))),kernel,'same'),'color',[0 .5 0])
hold on
plot(t,conv(full(sum(spiks(EI,:))),kernel,'same'))
axis tight
xlim(tLimPlot)


subaxis(4,1,3)
cla
% plot(f(f>=0),abs(spec(f>=0,:)))
% xlim([0 100])
imageNan(tt,f,squeeze(abs(spec)))
axis xy
% hold on
% boxCarLen=.3;
% h=plot(tt,conv(full(sum(spiks(EI,:))),ones(1,1e3*boxCarLen/dt),'same')/numE/boxCarLen,'color',[1 1 1]*.6,'linewidth',2);
% legend(h,'E firing rate')

% 
% subaxis(4,2,1,4)
% cla
% imagesc(reshape(output.I(1,[1:numE]+1),[sqrt(numE),sqrt(numE)]))


subaxis(4,1,1,4)
x=0:.5:50;
tSel=100:numel(t);
[Ne]=hist(sum(spiks(EI,tSel),2)/(numel(t(tSel))*dt*1e-3),x);
[Ni]=hist(sum(spiks(~EI,tSel),2)/(numel(t(tSel))*dt*1e-3),x);
% hold on
% bar(x,Ne/numE)
% bar(x,Ni/numI,'facecolor',[0 .5 0])
% xlim([3 max(x)])
% legend('E','I')
plot(x,[Ne(:)/numE,Ni(:)/numI])
xlabel('Firing rate')
xlim([3 max(x)])
legend('E','I')
xlabel('Firing rate')
title('histogram of F-rates')

fsz=14;
S_out=output.S;
S_orig=output.S_orig;
figure(11)
clf
set(11,'position',[150 150 1200 500])
subaxis(1,3,1)
imagesc(S_orig(1:numE,1:numE))
axis image
h=colorbar;
set(h,'location','southoutside')
title('original S','fontsize',fsz)
subaxis(1,3,2)
imagesc(S_out(1:numE,1:numE))
axis image
h=colorbar;
set(h,'location','southoutside')
title('new S','fontsize',fsz)
subaxis(1,3,3)
dS=S_out-S;
dS=dS(1:numE,1:numE);
imagesc(dS,maxabs(dS)*.1)
axis image
h=colorbar;
set(h,'location','southoutside')
title('\DeltaS','fontsize',fsz)

savefig=0;
if savefig
  fname='scholarp_multip_norm_10000A'
%   fname='Izh_STDP_hard_bound';
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs',[fname '_output']),'-png','-pdf','-transparent',10);
  export_fig(fullfile('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/output_figs',[fname '_Synap']),'-png','-pdf','-transparent',11);
end
