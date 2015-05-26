function [ connMat, connMat_sep, R_e, R_i ] = genConnections(numNeur, dR, sig, p, strength)
%[ connMat, connMat_sep, R_e, R_i ] = genConnections(numNeur, sig, p)
%   assumes 4:1 ratio of excitatory and inhibitory neurons
% 
%%%%INPUT%
%   numNeur:  the number of excitatory neurons along an edge of a square
%   lattice.
%   dR:       length of a lattice edge (in m)
%   sig:      standard deviations of (circular) gaussian connection
%   profiles (contains 4 values: EE, II, IE, EI)
%   p:        connection probability; fraction of possible sending neurons
%   that project to each target neuron (contains 4 values: EE, II, IE, EI)
%   strength: total conductance for each target neuron (contains 4 values: EE, II, IE, EI)
% 
%%%%OUTPUT%
%   connMat:  full connection matrix containing synaptic weigths
%   connMat_sep:  cell array containing the four sub matrices separatesly
%   R_e:      positions of excitatory neurons
%   R_i:      positions of inhibitory neurons

%% setting up neuron positions
% excitatory neurons
Rlin=0:dR:(numNeur-1)*dR;
% Rlin=Rlin-mean(Rlin);


[Rdum1,Rdum2]=meshgrid(Rlin);
R_e=nan(numel(Rdum1),2);
R_e(:,1)=Rdum1(:);
R_e(:,2)=Rdum2(:);

% inhibitory neurons
numNeurI=size(R_e,1)/4;
numNeurIlin=ceil(sqrt(numNeurI));

Rlin_i=linspace(0,numNeur*dR,numNeurIlin+1);
% Rlin_i=Rlin_i-mean(Rlin_i)
Rlin_i=Rlin_i(1:end-1);
[Rdum1,Rdum2]=meshgrid(Rlin_i);
R_i=nan(numel(Rdum1),2);
R_i(:,1)=Rdum1(:);
R_i(:,2)=Rdum2(:);


%% generating connection PDFs

EE=1;
II=2;
IE=3;
EI=4;

gridSize=numNeur*dR;

pdfConEE=connec_prob(sig(EE),gridSize,R_e);
pdfConII=connec_prob(sig(II),gridSize,R_i);
pdfConIE=connec_prob(sig(IE),gridSize,R_i,R_e);
pdfConEI=connec_prob(sig(EI),gridSize,R_e,R_i);

%% generate adjacency matrices

[A_EE]=sampleConn(pdfConEE,p(EE));
[A_II]=sampleConn(pdfConII,p(II));
[A_IE]=sampleConn(pdfConIE,p(IE));
[A_EI]=sampleConn(pdfConEI,p(EI));

%% transform to synaptic weights
connMat_sep=cell(2);
connMat_sep{EE}=bsxfun(@rdivide, A_EE, sum(A_EE,2))*strength(EE);
connMat_sep{II}=bsxfun(@rdivide, A_II, sum(A_II,2))*strength(II);
connMat_sep{IE}=bsxfun(@rdivide, A_IE, sum(A_IE,2))*strength(IE);
connMat_sep{EI}=bsxfun(@rdivide, A_EI, sum(A_EI,2))*strength(EI);

connMat=[connMat_sep{EE} connMat_sep{IE}; connMat_sep{EI} connMat_sep{II}];

connMat_sep=connMat_sep([EE, IE; EI, II]);

end

function [pdfCon]=connec_prob(sig, gridSize, R_send, R_receive)

if nargin<4
  R_receive=R_send;
  selfFlag=1;
else
  selfFlag=0;
end

x_send=R_send(:,1);
y_send=R_send(:,2);

x_receive=R_receive(:,1);
y_receive=R_receive(:,2);

DX=abs(bsxfun(@minus,x_receive,x_send'));
DY=abs(bsxfun(@minus,y_receive,y_send'));


if numel(gridSize)==1
  gridSize=[gridSize gridSize];
end
DX=min(DX,gridSize(1)-DX);
DY=min(DY,gridSize(2)-DY);


% gaussian connection pdf
pdfX=exp(-DX.^2/(2*sig^2));
pdfY=exp(-DY.^2/(2*sig^2));
pdfCon=pdfX.*pdfY;

if selfFlag
  pdfCon(logical(eye(size(pdfCon))))=0; % no self connections
end
pdfCon=bsxfun(@rdivide, pdfCon, sum(pdfCon,2));
end

function [connMat]=sampleConn(pdfCon,pCon)

[numReceive,numSend]=size(pdfCon);  
cpdf=cumsum(pdfCon,2);

%number of incoming connections per neuron (actually number of samples, therefore maximum)
numCon=round(pCon*numSend);

connMat=zeros(size(pdfCon));
dumsamp=rand(numCon, numReceive);
for n=1:numCon
  dum=bsxfun(@lt, dumsamp(n,:).', cpdf);
  [~, con]=max(dum,[],2);
  conidx=sub2ind([numReceive, numSend], [1:numReceive]', con);
  connMat(conidx)=connMat(conidx)+1;
end
end
