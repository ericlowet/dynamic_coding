function [ conMat, cfg, conMatSep] = genOrientColumns(cfg)
% [ connMat, connMat_sep, R_e, R_i ] = genOrientColumns(cfg)
% 
% 
% cfg contains the fields used:
% .gridSz:  number of hypercolumns along one axis of the square grid 
%           (total number of hyper columns will be the square of this number)
%           (default = 4)
% 
% .numOrient: number of orientation columns in every hyper column
%           (default = 4)      
% 
% .numE:    Number of excitatory neurons in every orientaiton column.
%           Number of inhibitory neurons will be 25% of this
%           (default = 40)
% 
% Connections WITHIN hypercolumns
% .numConIntra: [nEE nEI nIE nII] Number of incoming connections on every
%               neuron 
% 
% .strConIntra: [sEE sEI sIE sII] Strength of every incoming connection.
%               (default = [1e-2 5e-3 5e-2 5e-2])
% 
% .sig_theta:   standard deviation in radians of gaussian that determines
%               the probability for two neurons with preferred orientation
%               difference to be connected.
%               (note: also used for connections between hyper columns)
%               (default =  pi/8)
% 
% Connections BETWEEN hypercolumns
% .numConInter: [nEE nEI nIE nII] Number of incoming connections on every
%               neuron 
% 
% .strConInter: [sEE sEI sIE sII] Strength of every incoming connection.
%               (default = [1e-2 5e-3 5e-2 5e-2])
% 
% .sig_dist:    [dEE dEI dIE dII] standard deviation in grid units of
%               gaussian that determines probability of connections between
%               hypercolumns
%               (default = [1 1 .5 0])

%% parsing input
if isfield(cfg,'gridSz')
  gridSz=cfg.gridSz;
else
  gridSz=4;
  cfg.gridSz=gridSz;
end

if isfield(cfg,'numOrient')
  numOrient=cfg.numOrient;
else
  numOrient=4;
  cfg.numOrient=numOrient;
end

if isfield(cfg,'numE')
  numE=cfg.numE;
else
  numE=40;
  
  cfg.numE=numE;
end

if isfield(cfg,'numConIntra')
  numConIntra=cfg.numConIntra;
else
  numConIntra=[.5 .25 .1 .1].*numE;
  cfg.numConIntra=numConIntra;
end

if isfield(cfg,'strConIntra')
  strConIntra=cfg.strConIntra;
else
  strConIntra=[1e-2 5e-3 5e-2 5e-2];
  cfg.strConIntra=strConIntra;
end

if isfield(cfg,'numConInter')
  numConInter=cfg.numConInter;
else
  numConInter=[.5 .25 0 0].*numE*.5;
  cfg.numConInter=numConInter;
end

if isfield(cfg,'strConInter')
  strConInter=cfg.strConInter;
else
  strConInter=[1e-2 5e-3 5e-2 5e-2];
  cfg.strConInter=strConInter;
end

if isfield(cfg,'sig_theta')
  sig_theta=cfg.sig_theta;
else
  sig_theta=pi/8;
  cfg.sig_theta=sig_theta;
end

if isfield(cfg,'sig_dist')
  sig_dist=cfg.sig_dist;
else
  sig_dist=[1 1 .5 0];
  cfg.sig_dist=sig_dist;
end




%% setting up preferred orientations and spatial locations
numI=ceil(.25*numE);
x=[1:gridSz];%-(gridSz+1)/2
[X,Y]=meshgrid(x);
R=[X(:),Y(:)];

theta=linspace(0,pi,numOrient+1);
theta=theta(1:end-1);

conMatSep=cell(gridSz,gridSz);

numE_orient=numE/numOrient;
numI_orient=numI/numOrient;

numI_orient=ceil(numI_orient);
numI=numI_orient*numOrient;


%% connection WITHIN hyper column
% E -> E based on orientation
% E -> I randomly
% I -> E,I randomly
%
% numConIntra:  [nEE nEI nIE nII]






th_idx=repmat(1:numOrient,numE_orient,1);
th_idx=th_idx(:);

dist=0.5*angle(bsxfun(@times,exp(2i*theta),exp(2i*theta)'));
p_theta=exp(-0.5*dist.^2/sig_theta);
pdfCon_theta=p_theta(th_idx,th_idx);
pdfCon_theta=bsxfun(@rdivide,pdfCon_theta,sum(pdfCon_theta,2)); % normalize for /incoming/ connecitons
  

for hypColIdx=1:gridSz^2
  % EE
  conEE=strConIntra(1)*sampleCon(pdfCon_theta,numConIntra(1));
  
  % EI
  pdfCon=ones(numI,numE)/numE;
  conEI=strConIntra(2)*sampleCon(pdfCon,numConIntra(2));
  
  % IE
  pdfCon=ones(numE,numI)/numI;
  conIE=strConIntra(3)*sampleCon(pdfCon,numConIntra(3));
  
  % IE
  pdfCon=ones(numI,numI)/numI;
  conII=strConIntra(4)*sampleCon(pdfCon,numConIntra(4));
  
  conMatSep{hypColIdx,hypColIdx}=[conEE conIE; conEI conII];
end

%% Connections BETWEEN hypercolumns
% all connections based on both orientation and distance

gridLoc=meshgrid([1:gridSz]/(gridSz+1));
distx=(gridSz+1)*0.5/pi*angle(bsxfun(@times,exp(2i*gridLoc(:)*pi),exp(2i*gridLoc(:)*pi)'));
gridLoc=gridLoc.';
disty=(gridSz+1)*0.5/pi*angle(bsxfun(@times,exp(2i*gridLoc(:)*pi),exp(2i*gridLoc(:)*pi)'));
dist=sqrt(distx.^2+disty.^2');
nN=[numE, numI, numE, numI];
nNSnd=[numE, numE, numI, numI];


gridLocIdx=repmat(1:gridSz,gridSz,1);
gridLocIdx=gridLocIdx(:);
p_spat=cell(numel(sig_dist),1);
p_th=p_spat;
for n=1:numel(sig_dist)
  p_spat{n}=exp(-.5*dist.^2/sig_dist(n));
  p_spat{n}(logical(eye(gridSz)))=0;
  p_spat{n}(isnan(p_spat{n}(:)))=0;
  
  th_idx1=repmat(1:numOrient,nN(n)/numOrient,1);
  th_idx1=th_idx1(:);
  th_idx2=repmat(1:numOrient,nNSnd(n)/numOrient,1);
  th_idx2=th_idx2(:);
  p_th{n}=p_theta(th_idx1,th_idx2);
  p_th{n}=bsxfun(@rdivide,p_th{n},sum(p_th{n},2));
end




for hypColIdx1=1:gridSz^2
  numInterCon=cell(4,1);
  for n=1:numel(p_spat)    
    p_colmn=p_spat{n}((hypColIdx1),1:gridSz^2);
    if sum(p_colmn(1,:))>0
      p_colmn=p_colmn/sum(p_colmn(1,:));
      numInterCon{n}=full(sampleCon(p_colmn,numConInter(n)));
    else
      numInterCon{n}=zeros(1,gridSz^2);
    end    
  end  
  for hypColIdx2=1:gridSz^2
    if hypColIdx2~=hypColIdx1
      conDum=cell(4,1);
      for n=1:numel(numInterCon)
        conDum{n}=strConInter(n)*sampleCon(p_th{n},numInterCon{n}(hypColIdx2));
      end
      conMatSep{hypColIdx1,hypColIdx2}=[conDum{1} conDum{3}; conDum{2} conDum{4}];
    end
  end  
end

conMat=cell(gridSz^2,1);
for n=1:numel(conMat)
  conMat{n}=cat(2,conMatSep{n,:});
end
conMat=cat(1,conMat{:});
end

function [connMat]=sampleCon(pdfCon,nCon)

[nReceive,nSend]=size(pdfCon);
cpdf=cumsum(pdfCon,2);


connMat=sparse(zeros(size(pdfCon)));
dumsamp=rand(nCon, nReceive);
for n=1:nReceive
  [N,bin]=histc(dumsamp(:,n),[0 cpdf(n,:)]);
  connMat(n,bin)=N(bin);
end
end


