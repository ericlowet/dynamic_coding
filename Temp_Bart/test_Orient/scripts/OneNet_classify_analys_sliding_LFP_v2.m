basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/21_shiftMs';
basedir='/home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/classification/28';
addpath /home/mrphys/bargip/GIT/Dynamic_Coding/Temp_Bart/test_Orient/scripts
a=dir(fullfile(basedir,'N20*.mat'));
fnames={a.name};

subDir=[];

label='LFP';
label='Rate';
savDir=fullfile(basedir,['xCorrMat_' label '_v2']);
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

for pattIt=1:numel(fnames)
  disp(fnames{pattIt})
  load(fullfile(basedir,fnames{pattIt}))
  if pattIt==1
    dt=output.input.dt;    
    winLen=cfgms.ms_interval/dt;
  end
  
  rateDum=output.LFP;
  if strcmpi(label,'rate')
  rateDum=output.input.LFPkernel'*output.spiks;
  end
  rateDum=[rateDum(:,1)*0 rateDum];
  rateDum(isnan(rateDum))=0;
  rateDum=bg_reshape_overlap(rateDum, winLen, winLen, 2);
  
  rateDum=rateDum(:,:,srtIdx);
%   rateDum=rateDum(:,:,1:4:end);
  
  
  if pattIt==1
    rate=nan([size(rateDum) numel(fnames)]);
  end
  
  rate(:,:,:,pattIt)=rateDum;
end
rate=rate(:,:,:);

%%

winLen=50/dt;
winStarts=(1:winLen/2:(cfgms.ms_interval/dt-winLen+1));
numWin=numel(winStarts);

Xc=cell(1,numWin);
for winN=1:numWin
  rateDum=rate(:,winStarts(winN)-1+(1:winLen),:);
    tic
    [Xc{winN}] = slidXcorr(rateDum,1/dt,10/dt);
    toc
end

%%
numSac=cfgms.simulation_time./cfgms.ms_interval;

stimIdx=clp(repmat(phiIdx(srtIdx),1,numel(fnames)));
connIdx=clp(repmat(1:numel(fnames),numSac,1));

numStim=numel(unique(stimIdx));
numConn=numel(unique(connIdx));

combIdx=nan(size(connIdx));
count=0;
fac=1;
for n=1:numStim
  for m=1:numConn
    count=count+1;
    sel= stimIdx == n;
    sel= sel & connIdx ==m;
    combIdx( sel )=count;
  end
end
%%


for winN=1:numWin
  % decoding stimuli
  for cnI=1:numConn
    selTrials=connIdx==cnI;
    
    XcD=Xc{winN}(:,:,selTrials);
    XcD(isnan(XcD))=0;
    
    nT=size(XcD,3);
    nN=size(XcD,1);
    classVecS=cell(1,3);
    sel=repmat(logical(eye(nN)),1,1,nT);
    
    %fRates
    classVecS{1}=reshape(XcD(sel),nN,[]);
    
    %synchrony
    sel=repmat(triu(true(nN),1),1,1,nT);
    classVecS{2}=reshape(XcD(sel),nN*(nN-1)/2,[]);
    
    %lags
    XcD=permute(XcD,[2 1 3]);
    classVecS{3}=reshape(XcD(sel),nN*(nN-1)/2,[]);
    
    
    %   classVecS{3}=sign(classVecS{3});
    
    P=cell(1,1);
    P_ML=cell(1,1);
    tic
    [P{1},P_ML{1}]=xValPCAMnrClas(classVecS{1},stimIdx(selTrials),0,5);
    %   obj=fitcdiscr(classVecS{2}.',combIdx);
    [P{2},P_ML{2}]=xValPCAMnrClas(classVecS{2},stimIdx(selTrials),0,5);
    %   obj=fitcdiscr(classVecS{3}.',combIdx);
    [P{3},P_ML{3}]=xValPCAMnrClas((classVecS{3}),stimIdx(selTrials),0,5);
    [P{4},P_ML{4}]=xValPCAMnrClas([classVecS{1}; classVecS{3}],stimIdx(selTrials),0,5);
    %   [P{1},P_ML{1}]=xValPCAMnrClas(cat(1,classVecS{1:3}),combIdx,15,5);
    toc
    
    %   cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_1.mat'])},classVecS{1},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_2.mat'])},classVecS{2},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_3.mat'])},classVecS{3},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    save(fullfile(savDir,['P' num2str(winN,'%02d') '_conn' num2str(cnI) '_1msSm_10mSLag.mat']),'P_ML','P')
  end


  % decoding connectivity

  for stI=1:numStim
    selTrials=stimIdx==stI;
    
    XcD=Xc{winN}(:,:,selTrials);
    XcD(isnan(XcD))=0;
    
    nT=size(XcD,3);
    nN=size(XcD,1);
    classVecS=cell(1,3);
    sel=repmat(logical(eye(nN)),1,1,nT);
    
    %fRates
    classVecS{1}=reshape(XcD(sel),nN,[]);
    
    %synchrony
    sel=repmat(triu(true(nN),1),1,1,nT);
    classVecS{2}=reshape(XcD(sel),nN*(nN-1)/2,[]);
    
    %lags
    XcD=permute(XcD,[2 1 3]);
    classVecS{3}=reshape(XcD(sel),nN*(nN-1)/2,[]);
    
    
    %   classVecS{3}=sign(classVecS{3});
    
    P=cell(1,1);
    P_ML=cell(1,1);
    tic
    [P{1},P_ML{1}]=xValPCAMnrClas(classVecS{1},connIdx(selTrials),0,5);
    %   obj=fitcdiscr(classVecS{2}.',combIdx);
    [P{2},P_ML{2}]=xValPCAMnrClas(classVecS{2},connIdx(selTrials),0,5);
    %   obj=fitcdiscr(classVecS{3}.',combIdx);
    [P{3},P_ML{3}]=xValPCAMnrClas((classVecS{3}),connIdx(selTrials),0,5);
    [P{4},P_ML{4}]=xValPCAMnrClas([classVecS{1}; classVecS{3}],connIdx(selTrials),0,5);
    %   [P{1},P_ML{1}]=xValPCAMnrClas(cat(1,classVecS{1:3}),combIdx,15,5);
    toc
    
    %   cd ~/GIT/Dynamic_Coding/Temp_Bart/tq
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_1.mat'])},classVecS{1},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_2.mat'])},classVecS{2},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    %   qsubfeval('torqueWrap',{'xValPCAMnrClas',fullfile(savDir,['P' num2str(winN,'%02d') '_3.mat'])},classVecS{3},combIdx(:),10,5,'memreq',memrq,'timreq',timrq);
    save(fullfile(savDir,['P' num2str(winN,'%02d') '_stim' num2str(stI) '_1msSm_10mSLag.mat']),'P_ML','P')
  end
end

%%
classPerf=nan(numWin,4,2);
P_stim=cell(4,numWin);
P_conn=P_stim;

for winN=1:numWin
  P_stim(:,winN)={zeros(numStim)};
  P_conn(:,winN)={zeros(numConn)};
  % decoding stimuli
  for cnI=1:numConn
    load(fullfile(savDir,['P' num2str(winN,'%02d') '_conn' num2str(cnI) '_1msSm_10mSLag.mat']))
    for nn=1:numel(P_ML)
      P_stim{nn,winN}=P_stim{nn,winN}+P_ML{nn};
    end
  end
  % decoding connections
  for stI=1:numStim
    load(fullfile(savDir,['P' num2str(winN,'%02d') '_stim' num2str(stI) '_1msSm_10mSLag.mat']))
     for nn=1:numel(P_ML)
      P_conn{nn,winN}=P_conn{nn,winN}+P_ML{nn};
    end
  end
     
  for nn=1:numel(P_ML)

    P_stim{nn,winN}=P_stim{nn,winN}./sum(P_stim{nn,winN}(:));
    P_conn{nn,winN}=P_conn{nn,winN}./sum(P_conn{nn,winN}(:));
    
    classPerf(winN,nn,1)=sum(diag(P_stim{nn,winN}));
    classPerf(winN,nn,2)=sum(diag(P_conn{nn,winN}));
  end
end
%%
figure(11)
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
legend('energy','norm Corr','Lag','E+L','chance')
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
h=legend('energy','norm Corr','Lag','E+L','chance');
set(h,'location','northwest')
title('Connections')
xlim([0 max(winStarts/2+25)])
ylim([1/size(P_conn{1},1)*.8 max(clp(classPerf(:,:,2)))])

mtit(label,'fontsize',14,'yoff',.05)
%%
% addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,['classPerf_v2_10msLag_' label]),'-pdf','-png',11)

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
  load(fullfile(basedir,fnames{pattIt}))
  Ss{pattIt}=output.S;
end

figure(20)
clf

set(gcf,'position',[150 150 1200 800])
set(gcf,'color','w')

for pattIt=1:numel(fnames)
  subaxis(2,3,pattIt)
  imagesc(log(Ss{pattIt}))
  axis image
end

%%
addpath ~/test/Matlab/export_fig2014b/
export_fig(fullfile(basedir,'connProfiles'),'-pdf','-png',20)
