function [F, p]=multiGranger(data, maxLag, data_independent, indivFlag, savFile)
% [F, p,  B_full, B_res, X]=multiGranger(data, maxLag, data_independent)
% 
%%% INPUT %%%
% data=numChan x numTim x numTrl
% 
% maxLag = order of AR model
% 
% data_independent (optional) = more channels that are not targets, only
%                               senders (e.g. known input to the system)
%                               These are assumed to have an instantaneous
%                               effect.
% 
% indivFlag (optional) =  if data_independent is given determines whether
%                         data_independent needs to be regressed on all
%                         channels (0); or only on individual channels (1).
%                         Note: in the latter case, data_independent needs
%                         to have equal number of channels as data.
%                         (default = 1 when dimensions match, 0 otherwise)
% 
%%% OUTPUT %%%
% 
% F = ratio of variances of residuals (>1 means G-causal)
% 
% p = p-value of F in F-test for assessing significance of G-causality
% 
% B_full = Array containing regression coefficients for the full model
% 
% B_res = cell array containing restricted regression coeffient matrices
% 
% X = auto-regression matrix (containing the lagged signals)
% 
% Bart Gips November 2015

if nargin>4 && exist(savFile,'file')
  F=[];
  p=[];
  warning([savFile, ' already exists. Not performing analysis'])
  return
else
  indivFlag=false;
  [numChan, numTim, numTrl]=size(data);
  if nargin >2
    if ~isempty(data_independent)
      [nChi,nTi,nTri]=size(data_independent);
      if (nTi~=numTim) || (nTri ~= numTrl)
        error('number of timepoints or number of trials does not match for independent sources')
      end
      if nChi == numChan
        indivFlag=true;
      end
    end
  end 
  
  nT=numTim*numTrl;
  
  % building auto-regression matrix
  X=zeros([numTim, numTrl, maxLag, numChan]);
  % add independent sources if present
  if nargin>2
    if indivFlag
      for n=1:numChan
        dumy=data(n,:,:);
        dumx=data_independent(n,:,:);
        b=dumx(:)\dumy(:);
        dumy(:)=dumy(:)-dumx(:)*b;
        data(n,:,:)=dumy;
      end
    elseif ~isempty(data_independent)
      X=cat(2,X,permute(data_independent(:,:),[2 3 1]));
    end
  end
  
  
  for n=1:maxLag
    X(n+1:end,:,n,:)=permute(data(:,1:end-n,:),[2 3 1]);
  end
  X=reshape(X,[nT],maxLag,[]);
  
  Xdum=[ones(nT,1) X(:,:)]; %add constant
  
  Y=reshape(data,[numChan, nT]).';
  
  % calculate regression coefficients
  B_full=[Xdum(:,:) ]\Y;
  
  % calculate residuals
  res_tot=Y-[Xdum(:,:)]*B_full;
  var_full= var(res_tot).';
  
  %% now iteratively leave one channel out
  res=nan([size(res_tot) numChan]);
  
  if nargout >3 || nargin > 4
    B_res=cell(numChan,1);
  end
  
  for n=1:numChan
    sel=true(numChan,1);
    sel(n)=false; % leave current channel out
    Xdum=[ones(nT,1) reshape(X(:,:,sel),nT,[])]; %add constant
    
    b=[Xdum(:,:) ]\Y;
    
    if nargout >3 || nargin > 4
      B_res{n}=b;
    end
    res(:,:,n)=Y-[Xdum(:,:)]*b;
  end
  var_restr= squeeze(var(res));
  
  %% GC measure
  %  by way of ratio of variances (F-statistic) with corresponding statistical test
  F=bsxfun(@rdivide,var_restr,var_full);
  
  % one-tailed test; since var_restr can never be smaller than var_full
  p=1-fcdf(F,nT-1,nT-1);
  
  if nargin>4
    if ~exist(fileparts(savFile),'file')
      mkdir(fileparts(savFile));
    end
    %   save(savFile,'F', 'p', 'B_full', 'B_res', 'X')
    save(savFile,'F', 'p')
  end
end