basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/20';
basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/22_randCol';
basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/28';
addpath /home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/scripts
a=dir(fullfile(basedir,'N20*.mat'));
fnames={a.name};

subDir=[];

label='LFP2';
label='Rate2';
savDir=fullfile(basedir,['TFR_' label]);
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

foi=10:10:120;
foi=5:5:120;

for pattIt=1:numel(fnames)
  disp(fnames{pattIt})
  load(fullfile(basedir,fnames{pattIt}))
  if pattIt==1
    dt=output.input.dt;
    winLen=cfgms.ms_interval/dt;
    ddt=10e-3;
    toi=ddt/2:ddt:cfgms.simulation_time/1e3;
    
    wLFac=round(numel(output.t)/numel(toi));
  end
  
  rateDum=output.LFP;
  if strncmpi(label,'rate',4)
    rateDum=output.input.LFPkernel'*output.spiks;
  end
  rateDum=[rateDum(:,1)*0 rateDum ];
  rateDum(isnan(rateDum))=0;
%   rateDum=ft_specest_tfr(rateDum,[output.input.dt output.t]/1e3,'freqoi',foi,'timeoi',toi);
  rateDum=shiftdim(ft_specest_mtmconvol(rateDum,[output.input.dt output.t]/1e3,'freqoi',foi,'timeoi',toi,'timwin',100e-3*ones(size(foi)),'taper','hanning'),1);
  rateDum=bg_reshape_overlap(rateDum, winLen/wLFac, winLen/wLFac, 3);
  rateDum=permute(rateDum,[2 3 4 1]); % f x t x trial x chan
  rateDum=rateDum(:,:,srtIdx,:);
  
  cohDum=bsxfun(@times, conj(rateDum), permute(rateDum,[1 2 3 5 4]));
  sel=triu(true(8));
  cohDum=cohDum(:,:,:,sel(:));
  
  if pattIt==1
    coh=nan([size(cohDum) numel(fnames)]);
  end
  
  coh(:,:,:,:,pattIt)=cohDum;
end

save(fullfile(basedir,['TFRCoh_' label]),'coh','foi','toi');
%% generating indices for stimulus and connection pattern

numSac=cfgms.simulation_time./cfgms.ms_interval;

stimIdx=clp(repmat(phiIdx(srtIdx),1,numel(fnames)));
connIdx=clp(repmat(1:numel(fnames),numSac,1));


combIdx=nan(size(connIdx));
count=0;
fac=1;
for n=1:numel(unique(stimIdx))
  for m=1:numel(unique(connIdx))
    count=count+1;
    sel= stimIdx == n;
    sel= sel & connIdx ==m;
    combIdx( sel )=count;
  end
end
%%
load(fullfile(basedir,['TFRCoh_' label]));
% fsel=foi<=80;
% coh=coh(fsel,:,:,:,:);
% foi=foi(fsel);
% 
% fsel=4:4:32;
% coh=coh(fsel,:,:,:,:);
% foi=foi(fsel);

%%
numWin=size(coh,2);
numFreq=size(coh,1);

memrq=.5*1024^3;
timrq=60*60*5;

tmp=pwd;
if exist('/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/scripts/TFRClassHandle.mat','file')
  load /home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/scripts/TFRClassHandle.mat
else
  cd /home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/scripts
  batchid='TFRClass';
  compiledfun=qsubcompile(@dumFuncClass,'toolbox',{'stats'},'batchid',batchid);
  save TFRClassHandle compiledfun
end

cd(tmp)

for winN=1:numWin
  if ~exist(fullfile(savDir,['P' num2str(winN,'%02d') '_TFR_' label '.mat']),'file')
    cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
    qsubfeval(compiledfun,coh(:,winN,:,:,:),savDir,combIdx,winN,label,'memreq',memrq,'timreq',timrq);
  end
  %   dumFuncClass(coh(:,winN,:,:,:),savDir,combIdx,winN,label);
  
  
  %   P=cell(2,numFreq);
  %   P_ML=cell(2,numFreq);
  %   for freqN=1:numFreq
  %   disp(['freq ' num2str(freqN) ' out of ' num2str(numFreq)])
  %   XcD=permute(coh(freqN,winN,:,:,:),[4 3 5 1 2]);
  %   XcD=XcD(:,:);
  %   XcD(isnan(XcD))=0;
  %
  %   nT=size(XcD,2);
  %   nN=size(XcD,1);
  %   classVecS=cell(1,2);
  %
  %   %fRates
  %   sel=logical(eye(8));
  %   sel=sel(triu(true(size(sel))));
  %   sel=repmat(sel,1,nT);
  %   classVecS{1}=reshape(XcD(sel),[],nT);
  %
  %   %phase diffs
  %   sel=logical(triu(true(8),1));
  %   sel=sel(triu(true(size(sel))));
  %   sel=repmat(sel,1,nT);
  %   classVecS{2}=reshape(angle(XcD(sel)),[],nT);
  %
  %
  %
  %
  %   tic
  %   [P{1,freqN},P_ML{1,freqN}]=xValPCAMnrClas(classVecS{1},combIdx,0,5);
  %   [P{2,freqN},P_ML{2,freqN}]=xValPCAMnrClas(classVecS{2},combIdx,10,5);
  %   toc
  %   end
  %   save(fullfile(savDir,['P' num2str(winN,'%02d') '_TFR.mat']),'P_ML','P')
end

%%
numWin=size(coh,2);
numFreq=size(coh,1);
classPerf=nan(numFreq,numWin,2,2);
P_stim=cell(numFreq,numWin,2);
P_conn=P_stim;
Ps=zeros(numel(unique(stimIdx)));
Pc=zeros(numel(unique(connIdx)));
tic
for winN=1:numWin
  try
    load(fullfile(savDir,['P' num2str(winN,'%02d') '_TFR_' label '.mat']))
    
    
    for nn=1:size(P_ML,1);
      for ff=1:numFreq
        Pdum=P_ML{nn,ff};
        P_stim{ff,winN,nn}=Ps;
        P_conn{ff,winN,nn}=Pc;
        for n=1:max(combIdx)
          for m=1:max(combIdx)
            selx=combIdx==n;
            sely=combIdx==m;
            selCon=connIdx(selx);
            selStim=stimIdx(selx);
            conIx=selCon(1);
            stimIx=selStim(1);
            selCon=connIdx(sely);
            selStim=stimIdx(sely);
            conIy=selCon(1);
            stimIy=selStim(1);
            P_stim{ff,winN,nn}(stimIx,stimIy)=P_stim{ff,winN,nn}(stimIx,stimIy)+Pdum(n,m);
            P_conn{ff,winN,nn}(conIx,conIy)=P_conn{ff,winN,nn}(conIx,conIy)+Pdum(n,m);
          end
        end
        P_stim{ff,winN,nn}=P_stim{ff,winN,nn}./sum(P_stim{ff,winN,nn}(:));
        P_conn{ff,winN,nn}=P_conn{ff,winN,nn}./sum(P_conn{ff,winN,nn}(:));
        
        classPerf(ff,winN,nn,1)=sum(diag(P_stim{ff,winN,nn}));
        classPerf(ff,winN,nn,2)=sum(diag(P_conn{ff,winN,nn}));
      end
      
    end
  end
end
toc
%%
figure(sum(label))
clf
set(gcf,'position',[150 150 1200 800])
stimconLab={'Stimulus','Connection pattern'};
powerphaseLab={'Energy','Phase diff'};
repFac=1;
for stimconflag=1:2
  for powerphaseflag=1:2
    subaxis(2,2,stimconflag,powerphaseflag)
    cLims=[.2 nanmax(clp(classPerf(:,:,powerphaseflag,:)))];
    cLims=[.2 nanmax(clp(classPerf(:,:,:,:)))];
    imagesc(([1:size(classPerf(:,:,1),2)*repFac]-(repFac-1)*size(classPerf(:,:,1),2))*diff(toi(1:2)),foi,repmat(classPerf(:,:,powerphaseflag,stimconflag),[1 repFac]),cLims)
    colorbar
    if powerphaseflag==1
      title(stimconLab{stimconflag})
    else
      xlabel('Time (s)')
    end
    if stimconflag==1
      ylabel(sprintf([powerphaseLab{powerphaseflag} '\n Freq (Hz)']))
    end
    axis xy
  end
  
end
colormap hot
mtit(label(1:end),'fontsize',14,'yoff',.03)

% export_fig(fullfile(savDir,['ClassPerfTFR_' label]),'-pdf','-png',gcf)

%%
figure(12)
clf
set(gcf,'position',[150 150 1200 800])
subaxis(1,2,1)
plot(winStarts/2+25,classPerf(:,:,1))
hold on
plot([0 max(winStarts/2+25)],[1 1]/size(P_stim{1},1),'--k')
ax = gca;
ax.ColorOrderIndex = 1;
% plot(winStarts/2+25,classPerf(:,:,3),'--')
xlabel('window centre (ms)')
ylabel('class Perf')
legend('energy','norm Corr','Lag','chance')
title('Stimulus')
xlim([0 max(winStarts/2+25)])
ylim([1/size(P_stim{1},1)*.8 max(clp(classPerf(:,:,1)))])
subaxis(1,2,2)
plot(winStarts/2+25,classPerf(:,:,2))
hold on
plot([0 max(winStarts/2+25)],[1 1]/size(P_conn{1},1),'--k')
ax = gca;
ax.ColorOrderIndex = 1;
% plot(winStarts/2+25,classPerf(:,:,4),'--')

xlabel('window centre (ms)')
ylabel('class Perf')
h=legend(hlines,'energy','norm Corr','Lag','chance');
set(h,'location','northwest')
title('Connections')
xlim([0 max(winStarts/2+25)])
ylim([1/size(P_conn{1},1)*.8 max(clp(classPerf(:,:,2)))])

mtit(label,'fontsize',14,'yoff',.05)
%%
addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,['classPerf_10msLag_' label]),'-pdf','-png',12)

%% plot histograms
lab={'energy','correlation','lag'};
figure(31)
clf
set(gcf,'position',[100 150 1800 600])
for nn=1:size(P_stim,1)
  for kk=1:size(P_stim,2)
    subaxis(3,numWin,kk,nn, 'sv',.03,'sh',.02,'ml',.03,'mr',.01)
    imagesc(P_stim{nn,kk},minmax(cat(1,P_stim{:})))
    if nn==1;
      title([num2str((winStarts(kk)-1)/2+winLen/4) ' ms'])
    end
    if kk==1;
      ylabel(lab{nn})
    end
  end
  
end
colormap hot
mtit('stimulus','yoff',.05)

figure(32)
clf
set(gcf,'position',[100 150 1800 600])
for nn=1:size(P_conn,1)
  for kk=1:size(P_conn,2)
    subaxis(3,numWin,kk,nn, 'sv',.03,'sh',.02,'ml',.03,'mr',.01)
    imagesc(P_conn{nn,kk},minmax(cat(1,P_conn{:})))
    if nn==1;
      title([num2str((winStarts(kk)-1)/2+winLen/4) ' ms'])
    end
    if kk==1;
      ylabel(lab{nn})
    end
  end
  
end
colormap hot
mtit('connections','yoff',.05)

%%
addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,['classHist_' label]),'-pdf',31)
export_fig(fullfile(basedir,['classHist_' label]),'-pdf','-append',32)

%% plot connections

Ss=cell(1,numel(fnames));
for pattIt=1:numel(fnames)
  disp(fnames{pattIt})
  loaded=0;
  while ~loaded
  try
  load(fullfile(basedir,fnames{pattIt}))
  loaded =1;
  catch
    pause(1);
  end
  
  end
  Ss{pattIt}=output.S;
end

% adjacency matrices
totNE=cfg.numOrient*cfg.numE;
dumSel=false(totNE,1);
dumSel(1:cfg.numE)=true;
dumSel=repmat(dumSel,1,cfg.numOrient);
for n=2:cfg.numOrient
  dumSel(:,n)=circshift(dumSel(:,n-1),[cfg.numE 0]);
end

A=nan(cfg.numOrient,cfg.numOrient,numel(Ss));
k_in=nan(cfg.numOrient,numel(Ss));
k_out=k_in;
eigCentr=k_in;
for n=1:numel(Ss)
  dum=Ss{n}(1:totNE,1:totNE);
  A(:,:,n)=dumSel'*dum*dumSel;
  k_in(:,n)=sum(A(:,:,n),2);
  k_out(:,n)=sum(A(:,:,n),1);
  
  [V,D]=eig(A(:,:,n));
  eigCentr(:,n)=V(:,1);
  if all(V(:,1)<0)
    eigCentr(:,n)=-eigCentr(:,n);
  end
  
end
  
  


figure(20)
clf

set(gcf,'position',[150 150 1200 800])
set(gcf,'color','w')

for pattIt=1:numel(fnames)
  subaxis(2,3,pattIt,'mt',.05,'mb',.05,'ml',.05,'mr',.05)
  imagesc(log(Ss{pattIt}))
  axis image
  hcb=colorbar;
  if pattIt==numel(fnames)
    set(get(hcb,'ylabel'),'string','log(S)');
  end
end

figure(21)
clf

set(gcf,'position',[150 150 1200 800])
set(gcf,'color','w')

for pattIt=1:numel(fnames)
  subaxis(2,3,pattIt,'mt',.05,'mb',.05,'ml',.05,'mr',.05)
  imagesc(A(:,:,pattIt))
  axis image
  hcb=colorbar;
  if pattIt==numel(fnames)
    set(get(hcb,'ylabel'),'string','A');
  end
  title(['pattern ' num2str(pattIt)]);
end
mtit('Adjacency Matrices','yoff',0,'fontsize',12)

figure(22)
set(gcf,'position',[150 150 800 400])
set(gcf,'color','w')
clf
subaxis(1,2,1,'sh',.1)
plot(k_in,'linewidth',2)
axis tight
ylim(minmax(k_in)+[-1 1]*.1*max(k_in(:)))
ylabel('in-degree')
xlabel('orientation column #')

subaxis(1,2,2)
plot(eigCentr,'linewidth',2)
axis tight
ylim(minmax(eigCentr)+[-1 1]*.1*max(eigCentr(:)))
ylabel('Eigenvector Centrality')
xlabel('orientation column #')
legend(num2str([1:numel(fnames)]','Patt %d'))

% subaxis(1,3,2)
% plot(k_out)
% axis tight
% ylim(minmax(k_out)+[-1 1]*.1*max(k_out(:)))
% ylabel('out-degree')
% xlabel('orientation column #')

if ~exist(fullfile(basedir,'pattInfo'),'file')
  mkdir(fullfile(basedir,'pattInfo'))
end
export_fig(fullfile(basedir,'pattInfo','Adjacency_Matrix'),'-pdf',21)
export_fig(fullfile(basedir,'pattInfo','Graph_Theory'),'-pdf',22)

%%
addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,'connProfiles'),'-pdf','-png',20)
