clear all
addpath ~/MySims2/Dynamic_Coding/Temp_Bart/paruly
cd ~/MySims2/Dynamic_Coding/Temp_Bart/paruly
addpath('/home/common/matlab/fieldtrip/qsub')


% parameters RS
a=.02;   b=.2;
c=-65;    d=8;
parE=[a; b; c; d];
% parameters FS
a=.1;   b=.2;
c=-65;   d=2;
parI=[a; b; c; d];

%%%%%%%%%%%%%  TIME %%%%%%%%%%%%%%%
tLim=[0 2000]%60e3];
dt=.5;

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Connections %%%%%%
numNeurLin=31;
dR=50e-6;
% E->E, E->I, I->E, I->I
sig=[6    4  7 3]*dR;
%sig=[3 20 5 10]*dR;
%sig=[3 10 10 2]*dR;
p=[.2 .3 .3 .4];
strength=[.12 0.75 1 1];
[ connMat, connMat_sep, R_e, R_i ] = genConnPING(numNeurLin, dR, sig, p, strength);
numE=size(R_e,1);
numI=size(R_i,1);
S=connMat;
S_orig=S;

% figure,
% subplot(2,2,1)
% imagesc(S_orig(1: size(R_e,1)  , 1: size(R_e,1)  ))
% set(gca,'Clim', [0 0.005])
% 
% subplot(2,2,2)
% imagesc(S(size(R_e,1):  size(R_i,1)+ size(R_e,1), 1: size(R_e,1)  ))
% subplot(2,2,3)
% imagesc(S(1: size(R_e,1)  , size(R_e,1):  size(R_i,1)+ size(R_e,1)))
% subplot(2,2,4)
% imagesc(S(  size(R_e,1):  size(R_i,1)+ size(R_e,1)   , size(R_e,1):  size(R_i,1)+ size(R_e,1)))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%
% I=randn(numNeurLin,numNeurLin);
% sigI=4*dR;
% Rdum=min(R_e):dR:max(R_e);
% Rdum=Rdum-mean(Rdum);
% filtKern=exp(-Rdum.^2/(2*sigI^2));
% filtKern=bsxfun(@times,filtKern,filtKern');
% filtKern=filtKern/sum(filtKern(:));
% I=ifft2(fft2(I).*fft2(filtKern));
% I=I-min(I(:));
% I=I/mean(I(:))*4;
% I_orig=I;
%  figure,imagesc(I_orig)

net_size=numNeurLin;
 op=-pi:(2*pi)/(net_size-1):pi;
  op2=-pi:(2*pi)/(net_size-1):pi;
[X,Y] = meshgrid(op,op2); 
  conidx=sub2ind([net_size,net_size],round(net_size/2), round(net_size/2));
I = sqrt(circ_dist(X, X(conidx)).^2   + circ_dist(Y, Y(conidx)).^2 );
I=I - max(I(:))/2;% (I- median(I(:))  );
I= I./max(I(:));
I= (I./max(I(:))).*-0.5+4.7;
I_orig=I;
 figure,imagesc(I_orig)
set(gca,'Clim', [ 0 9])
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I=[I(:); zeros(size(R_i,1),1)]';
I=[tLim(1) I; tLim(2) I];

noise=I.*.04;
noise(:,1)=tLim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numNeur=size(R_e,1)+size(R_i,1);
V_init=repmat([c; b*c],1,numNeur);
V_init(1,:)=V_init(1,:)+randn(1,numNeur)*10;


numNeur=numE+numI;
EI=false(1,numNeur);
EI(1:numE)=true;

pars=[repmat(parE,1,numE) repmat(parI,1,numI)];
pars(3,1:numE)=pars(3,1:numE)+randn(1,numE)*3;
pars(4,1:numE)=pars(4,1:numE)+randn(1,numE)*1;
pars(1,numE+1:end)=pars(1,numE+1:end)+rand(1,numI)*.02;
pars(2,numE+1:end)=pars(2,numE+1:end)+rand(1,numI)*.02;


test_izhi='/home/coherence/erilow/MySims2/Temp_Bart/data/31net_circpos6.mat';
load(test_izhi)

%%%%%%%%%%%%
cfg=[];
cfg.tLim=tLim;
cfg.dt=dt;
cfg.a=pars(1,:);
cfg.b=pars(2,:);
cfg.c=pars(3,:);
cfg.d=pars(4,:);
cfg.S=output.S;%_orig;
cfg.S_orig =S_orig;
cfg.EI=EI;
cfg.STDP=0;
cfg.output='/home/coherence/erilow/MySims2/Temp_Bart/data/31net_circpos6_test1.mat';
cfg.saveInterval=500;%4e3;
cfg.verbose=0;
cfg.fullOutput=0;
cfg.A_STDP=[1 1]*1;



%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ~/MySims2/Temp_Bart/
cd ~/MySims2/Temp_Bart/batch

qsubfeval('Izh_network3',V_init,I,noise,cfg,'memreq',1024^3*9,'timreq',60^2*6);
%  [output,spikes]=Izh_network3(V_init,I,noise,cfg)
 % cfg.output=output;
cd '/home/coherence/erilow/MySims2/Temp_Bart/data/'





load(cfg.output)

figure,
imagesc(output.spiks)

rastersp= full(output.spiks);
[ signals ] = make_sig(rastersp(1:length(R_e),:),R_e,net_size);
xx=exp(i.*angle(signals(:,1:2:end)))';
 vv=(xx' *xx);
allcoh=single(angle(mean(vv,2)));
figure('COlor','w','Position', [ 300 300 300 200]),
imagesc(reshape(allcoh,31,31))
set(gca,'Clim', [ -0.5 0.5])

 amp=mean(abs(signals),2);
angs=(angle(signals));


figure('COlor','w','Position', [ 300 300 300 200]),
,imagesc(reshape(amp,net_size./1,net_size./1))
%figure,imagesc(I_orig)
figure,imagesc(output.S(1:length(R_e),1:length(R_e)))



