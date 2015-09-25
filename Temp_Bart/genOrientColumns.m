function [ conMat, idx, cfg, conMatSep] = genOrientColumns(cfg)
% [ conMat, idx, cfg, conMatSep] = genOrientColumns(cfg)
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
%           Number of inhibitory neurons will be 25% of this. Both will be
%           rounded up to match .numOrient
%           (default = 40)
% 
% Connections WITHIN hypercolumns
% .numConIntra: [nEE nEI nIE nII] Number of incoming connections on every
%               neuron 
% 
% .strConIntra: [sEE sEI sIE sII] Strength of every incoming connection.
%               (default = [5e-3 1e-2 4e-2 4e-2])
% 
% .sig_theta:   [E_intra I_intra E_inter I_inter] 
%               standard deviation in radians of gaussian that determines
%               the probability for two neurons with preferred orientation
%               difference to be connected. Separated for E-axons and
%               I-axons.
%               (note: also used for connections between hyper columns)
%               (default =  [pi/10 pi/5])
% 
% Connections BETWEEN hypercolumns
% .numConInter: [nEE nEI nIE nII] Number of incoming connections on every
%               neuron 
% 
% .strConInter: [sEE sEI sIE sII] Strength of every incoming connection.
%               (default = [5e-3 1e-2 4e-2 4e-2])
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

switch numel(gridSz)
  case 1
    gridSz=gridSz*[1 1];
  case 2
  otherwise
    error('gridSz should be a scalar or an array of 2 numbers')
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
end

if isfield(cfg,'numConIntra')
  numConIntra=cfg.numConIntra;
else
  numConIntra=ceil([.5 .5 .25 .25].*numE);
  cfg.numConIntra=numConIntra;
end

if isfield(cfg,'strConIntra')
  strConIntra=cfg.strConIntra;
else
%   [sEE sEI sIE sII]
  strConIntra=[5e-3 1e-2 4e-2 4e-2];
  cfg.strConIntra=strConIntra;
end

if isfield(cfg,'numConInter')
  numConInter=cfg.numConInter;
else
  numConInter=ceil([.5 .5 0 0].*numE*.5);
  cfg.numConInter=numConInter;
end

if isfield(cfg,'strConInter')
  strConInter=cfg.strConInter;
else
%   [sEE sEI sIE sII]
  strConInter=[5e-3 1e-2 4e-2 4e-2];
  cfg.strConInter=strConInter;
end

if isfield(cfg,'sig_theta')
  sig_theta=cfg.sig_theta;
  switch numel(sig_theta)
    case 1
      sig_theta=[2*pi 2*pi sig_theta 2*sig_theta];
    case 0
      sig_theta=[2*pi 2*pi pi/10 pi/5];
%     otherwise
%       error('sig_theta contains too many elements. Only two are required')
  end
else
  sig_theta=[2*pi 2*pi pi/10 pi/5];
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

cfg.numI=numI;

numHypColmn=prod(gridSz);

theta=linspace(0,pi,numOrient+1);
theta=theta(1:end-1);

conMatSep=cell(gridSz(1),gridSz(2));


nNeurRecv=[numE, numI, numE, numI];
nNeurSnd=[numE, numE, numI, numI];


%% connection probabilities

dist=0.5*angle(bsxfun(@times,exp(2i*theta),exp(2i*(theta))'));
p_theta_intra=cell(4,1);
p_theta_intra([1, 2])=deal({exp(-0.5*dist.^2/sig_theta(1))});
p_theta_intra([3, 4])=deal({exp(-0.5*dist.^2/sig_theta(2))});
p_theta_inter=cell(4,1);
p_theta_inter([1, 2])=deal({exp(-0.5*dist.^2/sig_theta(3))});
p_theta_inter([3, 4])=deal({exp(-0.5*dist.^2/sig_theta(4))});

gridLoc=meshgrid([1:gridSz(1)]/(gridSz(1)+1),[1:gridSz(2)]/(gridSz(2)+1));
X=gridLoc(:);
distX=(gridSz(1)+1)*0.5/pi*angle(bsxfun(@times,exp(2i*X*pi),exp(2i*X*pi)'));
gridLoc=gridLoc.';
Y=gridLoc(:);
distY=(gridSz(2)+1)*0.5/pi*angle(bsxfun(@times,exp(2i*Y*pi),exp(2i*Y*pi)'));
dist=sqrt(distX.^2+distY.^2');

p_spat=cell(numel(sig_dist),1);
[p_th_intra, p_th_inter]=deal(p_spat);
for n=1:numel(sig_dist)
  p_spat{n}=exp(-.5*dist.^2/sig_dist(n));
  p_spat{n}(logical(eye(prod(gridSz))))=0;
  p_spat{n}(isnan(p_spat{n}(:)))=0;
  
  th_idx1=repmat(1:numOrient,nNeurRecv(n),1);
  th_idx1=th_idx1(:);
  th_idx2=repmat(1:numOrient,nNeurSnd(n),1);
  th_idx2=th_idx2(:);
  p_th_intra{n}=p_theta_intra{n}(th_idx1,th_idx2);
  p_th_intra{n}=bsxfun(@rdivide,p_th_intra{n},sum(p_th_intra{n},2));
  p_th_inter{n}=p_theta_inter{n}(th_idx1,th_idx2);
  p_th_inter{n}=bsxfun(@rdivide,p_th_inter{n},sum(p_th_inter{n},2));
end

%% connection WITHIN hyper column
% all connections based on orientation; not location
%
% numConIntra:  [nEE nEI nIE nII]

for hypColIdx=1:numHypColmn
  conDum=cell(4,1);
  for n=1:numel(p_th_intra)
    conDum{n}=strConIntra(n)*sampleCon(p_th_intra{n},numConIntra(n));
  end
  
  conMatSep{hypColIdx,hypColIdx}=[conDum{1} conDum{3}; conDum{2} conDum{4}];
end

%% Connections BETWEEN hypercolumns
% all connections based on both orientation and distance


for hypColIdx1=1:numHypColmn
  numInterCon=cell(4,1);
  for n=1:numel(p_spat)    
    p_colmn=p_spat{n}((hypColIdx1),1:numHypColmn);
    if sum(p_colmn(1,:))>0
      p_colmn=p_colmn/sum(p_colmn(1,:));
      numInterCon{n}=full(sampleCon(p_colmn,numConInter(n)));
    else
      numInterCon{n}=zeros(1,numHypColmn);
    end    
  end  
  for hypColIdx2=1:numHypColmn
    if hypColIdx2~=hypColIdx1
      conDum=cell(4,1);
      for n=1:numel(numInterCon)
        conDum{n}=strConInter(n)*sampleCon(p_th_inter{n},numInterCon{n}(hypColIdx2));
      end
      conMatSep{hypColIdx1,hypColIdx2}=[conDum{1} conDum{3}; conDum{2} conDum{4}];
    end
  end
end

%%
conMat=cell(numHypColmn,1);
for n=1:numel(conMat)
  conMat{n}=cat(2,conMatSep{n,:});
end
conMat=cat(1,conMat{:});

%% output indices (i.e. neuron locations and orientation preference)
th_idxE=repmat(1:numOrient,numE,1);
th_idxE=th_idxE(:);
th_idxI=repmat(1:numOrient,numI,1);
th_idxI=th_idxI(:);

X=repmat(1:gridSz(1),(numE+numI)*numOrient*gridSz(2),1);
Y=repmat(1:gridSz(2),(numE+numI)*numOrient,gridSz(1));

EI=[ones(numE*numOrient,1); zeros(numI*numOrient,1)];

idx.X=X(:);
idx.Y=Y(:);
idx.theta=repmat([th_idxE; th_idxI],numHypColmn,1);
idx.EI=logical(repmat(EI,numHypColmn,1));


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


