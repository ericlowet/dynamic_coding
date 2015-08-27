cfg=[];
cfg.numOrient=6;
cfg.numE=60;
[ conMat, idx, cfg, conMatSep] = genOrientColumns(cfg);

%%
X=idx.X;
Y=idx.Y;
EI=idx.EI;
theta=idx.theta;

I=ones(size(theta))*.5;
I(~EI)=0;

sel=theta==2; 
selxy=X ==2 & Y==1;

I(sel & selxy)=I(sel & selxy)*2;
I(~sel & selxy)=I(~sel & selxy)*.5;

sel=theta==4; 
selxy=X ==3 & Y==3;
I(sel & selxy)=I(sel & selxy)*2;
I(~sel & selxy)=I(~sel & selxy)*4/5;

sel=theta==6; 
selxy=X ==1 & Y==4;
I(sel & selxy)=I(sel & selxy)*2;
I(~sel & selxy)=I(~sel & selxy)*4/5;

% sel=  X==2 | Y==4;
% 
% I(sel)=I(sel)*.5;
% 
sel= sqrt(((X-sum(minmax(X))/2).^2 + (Y-sum(minmax(Y))/2).^2))>1.6;

I(sel)=I(sel)*.5;

sel= sqrt(((X-sum(minmax(X))/2).^2 + (Y-sum(minmax(Y))/2).^2))>1;

I(sel)=I(sel)*.5;

I=I*5+(I>0)*5;

snr=2;

noise=(I>0).*(I(:))/snr;

figure(10)
set(gcf,'position',[200 200 800 600])
subaxis(2,2,1)
imagesc(reshape(mean(reshape(I(EI),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
colorbar
axis image
title('input strength')

subaxis(2,2,2,1)
imagesc(reshape(mean(reshape(noise(EI),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
colorbar
axis image
title('noise')

subaxis(2,2,1,2)
imagesc(.5*angle(reshape(mean(reshape(I(EI).*exp(1i*2*pi*theta(EI)/(cfg.numOrient)),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz)),[-pi pi]/2)
colorbar
axis image
title('input orientation')

subaxis(2,2,2,2)
imagesc(abs(reshape(mean(reshape(I(EI).*exp(1i*2*pi*theta(EI)/(cfg.numOrient)),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz)),[0 1])
colorbar
axis image
title('input orientation certainty (MVL)')


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



numNeur=numel(X);
V_init=repmat([c; b*c],1,numNeur);
V_init(1,:)=V_init(1,:)+randn(1,numNeur)*1;

pars=nan(4,numNeur);
for n=1:4
  pars(n,EI)=parE(n)*(1+.01*randn(1,sum(EI)));
  pars(n,~EI)=parI(n)*(1+.01*randn(1,sum(~EI)));
end


tLim=[0 2e3];
dt=.5;

input=[];
input.I=[tLim(:) [I(:)'; I(:)']];
input.noise = [ tLim(:) [noise(:)' ; noise(:)']];

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

% input.output=['~/GIT/Dynamic_Coding/Temp_Bart/explicit_Syn/N31_clock_snr2_2attrac.mat'];

[output,spikes]=Izh_network_TAH(input);

%



%%
plotLim=tLim;%[0 2]*1e3;

dt=output.input.dt;
EI=output.EI;

input=output.input;


ksig=2/dt;
tt=-5*ksig:5*ksig;
sphistkernel=exp(-tt.^2/(2*ksig^2));
sphistkernel=sphistkernel/sum(sphistkernel);
figure(2)
set(2,'position',[100 100 1200 800]) 
clf
subaxis(5,1,1)
hold all
plot(output.t,conv(full(sum(output.spiks(EI,:))),sphistkernel(:),'same'))
plot(output.t,conv(full(sum(output.spiks(~EI,:))),sphistkernel(:),'same'))
xlabel('time (ms)')
xlim([output.input.tLim(1) max(output.t)])
xlim(plotLim)
subaxis(5,1,2)
ns_plotspikes(spikes,gca, [], plotLim)
xlim([output.input.tLim(1) max(output.t)])
xlim(plotLim)


sphist=cat(1,full(sum(output.spiks(EI,:))),full(sum(output.spiks(~EI,:)))).';
[sp,f]= bg_fftwelch(sphist,1e3/dt,.5,1);
subaxis(5,1,3)
plot(f(f>=0),mean(abs(sp(f>=0,:,:)).^2,3))
xlim([0 100])
ylabel('power')
xlabel('freq (Hz)')

subaxis(5,4,1,4,1,2)
imagesc(reshape(mean(reshape(sum(output.spiks(EI,:),2)/((diff(output.I([1 end],1))+output.input.dt)/1e3),cfg.numE,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
colorbar
axis image
title('rate E')

subaxis(5,4,2,4,1,2)
imagesc(reshape(mean(reshape(sum(output.spiks(~EI,:),2)/((diff(output.I([1 end],1))+output.input.dt)/1e3),cfg.numI,cfg.gridSz^2)),cfg.gridSz,cfg.gridSz))
colorbar
axis image
title('rate I')

dum=nan(cfg.gridSz,cfg.gridSz,numel(unique(idx.theta)));
dumAng=dum;
for n=1:numel(unique(idx.theta))
  dum(:,:,n)=reshape(mean(reshape(sum(output.spiks(EI & idx.theta==n,:),2),[],cfg.gridSz^2)),cfg.gridSz,cfg.gridSz);
  dumAng(:,:,n)=exp(1i*2*pi*(n/(numel(unique(idx.theta)))));
end
th_dec=sum(dum.*dumAng,3);

subaxis(5,4,3,4,1,2)
imagesc(abs(th_dec))
axis image
colorbar

subaxis(5,4,4,4,1,2)
imagesc(angle(th_dec)/2,[-pi pi]/2)
axis image
colorbar
% colormap jet

figure(21)
imagesc(angle(th_dec)/2,[-pi pi]/2)
axis image
colormap hsv
