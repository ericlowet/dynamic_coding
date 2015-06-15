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
tLim=[0 6000]%60e3];
dt=.5;

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Connections %%%%%%
numNeurLin=31;
dR=50e-6;
% E->E, E->I, I->E, I->I
sig=[8    4  7 3]*dR;
%sig=[3 20 5 10]*dR;
p=[.4 .3 .3 .4];
strength=[.2 0.75 1 1];
[ connMat, connMat_sep, R_e, R_i ] = genConnPING(numNeurLin, dR, sig, p, strength);
numE=size(R_e,1);
numI=size(R_i,1);
S=connMat;
S_orig=S;
figure,
subplot(2,2,1)
imagesc(S_orig(1: size(R_e,1)  , 1: size(R_e,1)  ))
subplot(2,2,2)
imagesc(S(size(R_e,1):  size(R_i,1)+ size(R_e,1), 1: size(R_e,1)  ))
subplot(2,2,3)
imagesc(S(1: size(R_e,1)  , size(R_e,1):  size(R_i,1)+ size(R_e,1)))
subplot(2,2,4)
imagesc(S(  size(R_e,1):  size(R_i,1)+ size(R_e,1)   , size(R_e,1):  size(R_i,1)+ size(R_e,1)))




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
I= (I./max(I(:))).*2.5+5.4;
I_orig=I;
 figure('COlor','w'),
 subplot(1,1,1,'Fontsize',17)
 imagesc(I_orig)
colorbar
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I=[I(:); zeros(size(R_i,1),1)]';
I=[tLim(1) I; tLim(2) I];

noise=I.*.03;
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
pars(1,numE+1:end)=pars(1,numE+1:end)+rand(1,numI)*.025;
pars(2,numE+1:end)=pars(2,numE+1:end)+rand(1,numI)*.025;




%%%%%%%%%%%%
cfg=[];
cfg.tLim=tLim;
cfg.dt=dt;
cfg.a=pars(1,:);
cfg.b=pars(2,:);
cfg.c=pars(3,:);
cfg.d=pars(4,:);
cfg.S=S;
cfg.S_orig =S_orig;
cfg.EI=EI;
cfg.STDP=1;
cfg.output='/home/coherence/erilow/MySims2/Temp_Bart/data/31net_circneg6.mat';
cfg.saveInterval=1000;%4e3;
cfg.verbose=0;
cfg.fullOutput=0;
cfg.A_STDP=[0.00004 0.00004];



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
% figure,
% imagesc(output.spiks)
rastersp= full(output.spiks);
rastersp((rastersp==0)) =NaN;
figure('Color','w')
subplot(1,1,1,'Fontsize',17);
for ind=1:length(R_e)
hold on, plot( (1:length(rastersp(1,:)))./2000 ,  rastersp(ind,:)+ind, '.', 'Markersize',10)
end
for ind=length(R_e)+1 : size(rastersp,1)
hold on, plot( (1:length(rastersp(1,:)))./2000 ,  rastersp(ind,:)+ind, '.', 'Markersize',10,'COlor', [ 0 0.5 0])
end

S=output.S;

figure('Color','w')
subplot(1,2,2)
imagesc(S(1: size(R_e,1)  , 1: size(R_e,1)  ))
%set(gca,'Clim', [0.0005 0.004])
subplot(1,2,1)
imagesc(S_orig(1: size(R_e,1)  , 1: size(R_e,1)  ))
%set(gca,'Clim', [0.0005 0.004])

figure('Color','w')
subplot(1,1,1,'Fontsize',17);
imagesc(smooth2(S(1: size(R_e,1)  , 1: size(R_e,1)  ) - S_orig(1: size(R_e,1)  , 1: size(R_e,1)  ),2,2))
set(gca,'Clim', [-0.001 0.001])

%figure,plot(output.deltaS(:,1).*-1)
%hold on, plot(output.deltaS(:,2).*1,'r')




% subplot(2,2,2)
% imagesc(S(size(R_e,1):  size(R_i,1)+ size(R_e,1), 1: size(R_e,1)  ))
% subplot(2,2,3)
% imagesc(S(1: size(R_e,1)  , size(R_e,1):  size(R_i,1)+ size(R_e,1)))
% subplot(2,2,4)
% imagesc(S(  size(R_e,1):  size(R_i,1)+ size(R_e,1)   , size(R_e,1):  size(R_i,1)+ size(R_e,1)))






