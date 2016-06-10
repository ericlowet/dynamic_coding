addpath /home/mrphys/bargip/test/Matlab/te_matlab


basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/double5';


a=dir(fullfile(basedir,'*.mat'));
fnames={a.name};

subDir=[];

savDir=fullfile(basedir,'xCorrMat');
if ~exist(savDir,'file')
  mkdir(savDir)
  %   mkdir(fullfile(savDir,subDir))
end
savFname=fnames;
for n=1:numel(fnames)
  savFname{n}=savFname{n}(1:end-4); %cut off .mat
end



load(fullfile(basedir,'stim','decInp.mat'));
load(fullfile(basedir,'stim','aux.mat'));
[~,srtIdx]=sort(phiIdx);
if exist('idx','var')
  X=idx.X;
  Y=idx.Y;
  EI=idx.EI;
  theta=idx.theta;
  numOrient=numel(unique(theta));
end
%%
fnames=fnames(1);

for pattIt=1:numel(fnames)
  disp(fnames{pattIt})
  load(fullfile(basedir,fnames{pattIt}))
  if pattIt==1
    dt=output.input.dt;
    winLen=cfgms.ms_interval/dt;
  end
  
  rateDum=[output.spiks(EI,:) output.spiks(EI,1)*0];
  rateDum(isnan(rateDum))=0;
  rateDum=bg_reshape_overlap(rateDum, winLen, winLen, 2);
%   rateDum=rateDum(:,:,srtIdx);
  if pattIt==1
    rate=nan([size(rateDum) numel(fnames)]);
  end
  
  rate(:,:,:,pattIt)=rateDum;
end

winLen=100/dt;
winStarts=1:winLen/2:(cfgms.ms_interval/dt-winLen+1);
maxLag=winLen/2;

numWin=numel(winStarts);
%%
% P=cell(numWin,2,3);
timrq=120*60;
memrq=1*1024^3;
for winN=1%:numWin
  rateDum=rate(:,winStarts(winN)-1+(1:winLen),:);
  
  % zero-padding
  rateDum=[zeros([size(rateDum,1) maxLag size(rateDum,3)]) rateDum];
  rateDum=sparse(rateDum(:,:));
  
  asdf=SparseToASDF(rateDum,dt);
%   tic
%   [TE, TE1_lag] = ASDFTE(asdf,1:maxLag,4,4);
%   toc
  
  cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
  qsubfeval('torqueWrap',{'ASDFTE',fullfile(basedir,'TE_ho5',[num2str(pattIt,['P_' num2str(winLen*dt) 'ms_%d']) '_' num2str(winN,'%02d')])},asdf,1:maxLag,5,5,'memreq',memrq,'timreq',timrq);
  
end
%%
nN=160;

S=output.S(EI,EI);
S_intra=[S(1:nN,1:nN); S(nN+[1:nN],nN+[1:nN])];
S_inter=[S(nN+[1:nN],1:nN); S([1:nN],nN+[1:nN])];

kernSz=1;


winLen1=50/dt;
winStarts1=1:winLen1/2:(cfgms.ms_interval/dt-winLen1+1);

numWin=numel(winStarts1);
c=nan(numWin,2);
MI=nan(numWin,2);
numBinsHist=24;
for winN=1:numWin
  try
  load(fullfile(basedir,'TE_ho5',[num2str(pattIt,'P_50ms_%d') '_' num2str(winN,'%02d')]))
%   load(fullfile(basedir,'TE_ho5',[num2str(pattIt,'P_100ms_%d') '_' num2str(winN,'%02d')]))
  TE=outp{1};
  TE_intra=[TE(1:nN,1:nN); TE(nN+[1:nN],nN+[1:nN])];
  TE_inter=[TE(1:nN,nN+[1:nN]); TE([1:nN],nN+[1:nN])];
  
%   figure(2)
%   subaxis(3,5,winN)
%   imagesc(TE)
  TE(isnan(TE))=0;
  c(winN,1)=corr(clp(convn(TE_intra,ones(kernSz))),clp(convn(full(S_intra),ones(kernSz))),'type','spearman');
  c(winN,2)=corr(clp(convn(TE_inter,ones(kernSz))),clp(convn(full(S_inter),ones(kernSz))),'type','spearman');
  
  MI(winN,1)=miNorm(S_intra(:),TE_intra(:), {linspace(0,max(S_intra(:)),numBinsHist), linspace(min(TE_intra(:)),max(TE_intra(:)),numBinsHist)},0);
  MI(winN,2)=miNorm(S_inter(:),TE_inter(:), {linspace(0,max(S_inter(:)),numBinsHist), linspace(min(TE_inter(:)),max(TE_inter(:)),numBinsHist)},0);
  
  end
end


winLen2=100/dt;
winStarts2=1:winLen2/2:(cfgms.ms_interval/dt-winLen2+1);
numWin=numel(winStarts2);
c2=nan(numWin,2);
MI2=nan(numWin,2);
for winN=1:numWin
  try
%   load(fullfile(basedir,'TE_ho5',[num2str(pattIt,'P%d') '_' num2str(winN,'%02d')]))
  load(fullfile(basedir,'TE_ho5',[num2str(pattIt,'P_100ms_%d') '_' num2str(winN,'%02d')]))
  TE=outp{1};
  TE_intra=[TE(1:nN,1:nN); TE(nN+[1:nN],nN+[1:nN])];
  TE_inter=[TE(1:nN,nN+[1:nN]); TE([1:nN],nN+[1:nN])];
%   
%   figure(2)
%   subaxis(3,5,winN)
%   imagesc(TE)
  TE(isnan(TE))=0;
  c2(winN,1)=corr(clp(convn(TE_intra,ones(kernSz))),clp(convn(full(S_intra),ones(kernSz))),'type','spearman');
  c2(winN,2)=corr(clp(convn(TE_inter,ones(kernSz))),clp(convn(full(S_inter),ones(kernSz))),'type','spearman');
  
  MI2(winN,1)=miNorm(S_intra(:),TE_intra(:), {linspace(0,max(S_intra(:)),numBinsHist), linspace(min(TE_intra(:)),max(TE_intra(:)),numBinsHist)},0);
  MI2(winN,2)=miNorm(S_inter(:),TE_inter(:), {linspace(0,max(S_inter(:)),numBinsHist), linspace(min(TE_inter(:)),max(TE_inter(:)),numBinsHist)},0);
 
  end
end
%%
fsz =16;
figure(8)
clf
set(gcf,'position',[150 150 700 450])
subaxis(1,2,1,'mr',.02,'sh',.1,'mb',.1)
set(gca,'fontsize',fsz)
imagesc(log(output.S))
cb=colorbar('location','southoutside');
set(get(cb,'xlabel'),'string','log(S)')
axis image
colormap hot
% subaxis(2,2,2,1)
subaxis(1,2,2)
hold all
ax = gca;
set(ax,'ColorOrder', [0.8500    0.3250    0.0980; 0    0.4470    0.7410]);

plot(winStarts1*dt+winLen1/2*dt,c)

set(gca,'fontsize',fsz)

plot(winStarts2*dt+winLen2/2*dt,c2,'--')
ylim([0 1])
xlabel('window centre (ms)')
ylabel('correlation (TE\astS)')
hh=legend('Intra (50ms window)','Inter (50ms window)','Intra (100ms window)','Inter (100ms window)');
set(hh,'location','northwest')
% 
% subaxis(2,2,2,2)
% plot(winStarts1*dt+winLen1/2*dt,MI)
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(winStarts2*dt+winLen2/2*dt,MI2,'--')
% xlabel('window centre (ms)')
% ylabel('I(TE,S)')
% hh=legend('Intra (50ms window)','Inter (50ms window)','Intra (100ms window)','Inter (100ms window)');
% set(hh,'location','best')


%%
addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,'TE_summ'),'-transparent','-pdf','-eps','-png',8)