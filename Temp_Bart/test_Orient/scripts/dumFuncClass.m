function dumFuncClass(coh,savDir,combIdx,winN,label)

numFreq=size(coh,1);
P=cell(2,numFreq);
P_ML=cell(2,numFreq);
for freqN=1:numFreq
  disp(['freq ' num2str(freqN) ' out of ' num2str(numFreq)])
  XcD=permute(coh(freqN,1,:,:,:),[4 3 5 1 2]);
  XcD=XcD(:,:);
  XcD(isnan(XcD))=0;
  
  nT=size(XcD,2);
  nN=size(XcD,1);
  classVecS=cell(1,2);
  
  %fRates
  sel=logical(eye(8));
  sel=sel(triu(true(size(sel))));
  sel=repmat(sel,1,nT);
  classVecS{1}=reshape(XcD(sel),[],nT);
  
  %phase diffs
  sel=logical(triu(true(8),1));
  sel=sel(triu(true(size(sel))));
  sel=repmat(sel,1,nT);
  classVecS{2}=reshape(angle(XcD(sel)),[],nT);
  
  
  
  
  tic
  [P{1,freqN},P_ML{1,freqN}]=xValPCAMnrClas(classVecS{1},combIdx,0,5);
  [P{2,freqN},P_ML{2,freqN}]=xValPCAMnrClas(classVecS{2},combIdx,10,5);
  toc
end
save(fullfile(savDir,['P' num2str(winN,'%02d') '_TFR_' label '.mat']),'P_ML','P')
end