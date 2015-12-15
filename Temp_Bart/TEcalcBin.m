function [TE, TE_Norm]=TEcalcBin(spiks, lags)
% [TE_delay, TE_delayNorm]=TEcalcBin(spiks, lags)
%
% based on:
% Ito, Shinya, et al.
% "Extending transfer entropy improves identification of effective connectivity in a spiking cortical network model."
% PloS one 6.11 (2011): e27431.
%
% Essentially a replacement/extension of their transent.c function
% (I am not proficient with C, so have ported it to MATLAB)
% 
% In short: it assumes spiking is sparse, therefore it only visits all the 
% time bins worth visiting (i.e. with spikes) and then fills in the silent 
% time bins automatically.
%
% Bart Gips
% December 2015

numNeur=spiks{end}(1);
numBins=spiks{end}(2);

if nargin<2
  lags=1;
end
nLags=numel(lags);

TE=nan(numNeur,numNeur,nLags);
TE_Norm=TE;

for nN1=1:numNeur % target neuron nN1
  for nN2=[1:nN1-1 nN1+1:numNeur] % sending neuron nN2
    
    %keep histogram count [n1(t), n1(t-1), n2(t-d), 1:d]
    jointHist=zeros(2,2,2,nLags);
    % counting all possible patterns
    
    
    % 1) when target is spiking check if target spiked in t-1 or sender
    % spiked in t-d
    for nS=1:numel(spiks{nN1})
      % target spikes at t
      t_sp=spiks{nN1}(nS);
      lIdx=t_sp-lags;
      
      sendSpiksDum=spiks{nN2};
      sendSpiksDum=sendSpiksDum(sendSpiksDum<(t_sp+1) & sendSpiksDum>(lIdx(end)));
      
      if nS>1 % check if target neuron fired in t-1
        sp1=(t_sp-1)==spiks{nN1}(nS-1);
      else
        sp1=false;
      end
      
      sp2=any(bsxfun(@eq,lIdx(:),sendSpiksDum),2);
      
      for lag=1:nLags
        if t_sp>lags(lag) % only if we know for sure that there was or was no spike (i.e. we have data)
          jointHist(2,sp1+1,sp2(lag)+1,lag)=jointHist(2,sp1+1,sp2(lag)+1,lag)+1;
        end
      end
      
      % 2) move over 1 timestep, i.e. target now spikes at t-1;
      % Secondly next bin should be empty. Otherwise it
      % will be covered by the step 1 (spike at t)
      % (only when next bin is in the data.)
      if spiks{nN1}(nS)+1 < numBins && (nS==numel(spiks{nN1}) || spiks{nN1}(nS+1) ~= spiks{nN1}(nS)+1)
        
        
        sp2=any(bsxfun(@eq,lIdx(:)+1,sendSpiksDum),2);
        for lag=1:nLags
          if t_sp>lags(lag) % only if we know for sure that there was or was no spike (i.e. we have data)
            jointHist(1,2,sp2(lag)+1,lag)=jointHist(1,2,sp2(lag)+1,lag)+1;
          end
        end
      end
      
    end
    
    % 3) when sender is spiking at t-d and target is silent in both t and
    % t-1 bins [those two options are covered in 1) and 2) ]
    for nS=1:numel(spiks{nN2})
      t_sp=spiks{nN2}(nS)+lags;
      
      targSpiksDum=spiks{nN1};
      targSpiksDum=targSpiksDum(targSpiksDum>=t_sp(1) & targSpiksDum<t_sp(end));
      
      sp_sel=~any(bsxfun(@eq,t_sp,targSpiksDum(:))) & ~any(bsxfun(@eq,t_sp-1,targSpiksDum(:))) & t_sp<numBins;
      jointHist(1,1,1,~sp_sel)=jointHist(1,1,1,~sp_sel)+1;
      jointHist(1,1,2,sp_sel)=jointHist(1,1,2,sp_sel)+1;
    end
    % 4) all other time bins belong to no spikes at all
    jointHist(1,1,1,:)=bsxfun(@minus,numBins,lags*0).'-permute(sum(sum(sum(jointHist,1),2),3),[4 1 2 3]);
    % note: normally the 0 in lags*0 should not be there. But because of
    % sparseness of firing, changing it to lags*1 will lead to a slight systematic
    % increase in the TE with increased lag. In cases of low SNR this
    % amounts to the location of maximum TE being biased towards high lags.
    % Changing this to lags*0 essentially sets every pattern that depends
    % on spikes outside of the data to the non-informative "0" pattern
    % (i.e. total silence). This is similar (but not equal) to 0-padding
    % the data up to maxLags.
    
    jointProb=bsxfun(@rdivide,jointHist,sum(sum(sum(jointHist,1),2),3));
    condProb1=sum(jointProb,3);
    condProb1=bsxfun(@rdivide,condProb1,sum(condProb1));
    condProb2=bsxfun(@rdivide,jointProb,sum(jointProb));
    TE_delay_log=bsxfun(@minus,log2(condProb2),log2(condProb1));
    TE_delaydum=bsxfun(@times,TE_delay_log,jointProb);
    TE_delaydum(~jointProb)=0; % 0 * log(x) = 0; also if x=inf
    TE(nN1,nN2,:)=permute(sum(sum(sum(TE_delaydum,1),2),3),[4 1 2 3]);
    
    % Since TE(x->y)=H(y|y-)-H(y|y-,x-)  [where x- means past of x]
    % normalize by conditional entropy: H(y|y-). 
    
    % same for all delays if "lags*0" is used above
    pyx=sum(jointProb(:,:,:,:),3);
    pyGivx=bsxfun(@rdivide,pyx,sum(pyx,1));
    HyGivx=-bsxfun(@times,pyx,log2(pyGivx));
    HyGivx=nansum(nansum(HyGivx,1),2);
    
    TE_Norm(nN1,nN2,:)=TE(nN1,nN2,:)./permute(HyGivx,[1 2 4 3]);
    
  end
end
end
