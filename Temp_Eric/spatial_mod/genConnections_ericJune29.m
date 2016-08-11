function [ connMat, connMat_sep, R_e, R_i ] = genConnections_ericV2(numNeur, dR, sig, p, strength,E2num,I2num)
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
%%%%OUTPUT%
%   connMat:  full connection matrix containing synaptic weigths
%   connMat_sep:  cell array containing the four sub matrices separatesly
%   R_e:      positions of excitatory neurons
%   R_i:      positions of inhibitory neurons
%   R_e2:      positions of excitatory neurons
%   R_i2:      positions of inhibitory neurons
%% setting up neuron positions

%% excitatory neurons
Rlin=0:dR:(numNeur-1)*dR;
numNeurE = numNeur*numNeur;
[Rdum1,Rdum2]=meshgrid(Rlin);
R_e=nan(numel(Rdum1),2);
R_e(:,1)=Rdum1(:);
R_e(:,2)=Rdum2(:);

% V2 E neurons
numNeurE2=E2num;
midE=round(size(Rdum1,1)/2);
R_e2(:,1)=ones(numNeurE2,1).*Rdum1(midE,midE) ;R_e2(:,2)=ones(numNeurE2,1).*Rdum2(midE,midE) ;

%% inhibitory neurons
numNeurI=size(R_e,1)/4;
numNeurIlin=ceil(sqrt(numNeurI));
Rlin_i=linspace(0,numNeur*dR,numNeurIlin+1);
Rlin_i=Rlin_i(1:end-1);
[Rdum1,Rdum2]=meshgrid(Rlin_i);
R_i=nan(numel(Rdum1),2);
R_i(:,1)=Rdum1(:);R_i(:,2)=Rdum2(:);
% V2 I neurons
numNeurI2=I2num;
R_i2(:,1)=ones(numNeurI2,1).*Rdum1(midE,midE) ;R_i2(:,2)=ones(numNeurI2,1).*Rdum2(midE,midE) ;

totalN=numNeurE+numNeurE2+numNeurI+numNeurI2;
E_ind=1:numNeurE;
E2_ind=numNeurE+1:numNeurE+numNeurE2;
I_ind=(E2_ind(end)+1):(E2_ind(end)+numNeurI);
I2_ind=(E2_ind(end)+numNeurI)+1:totalN;

%% generating connection PDFs
EE=1;II=2;IE=3;EI=4;EE2=5;
E2E2=6;E2I2=7;I2I2=8;I2E2=9;

gridSize=numNeur*dR;
pdfConEE=connec_prob(sig(EE),gridSize,R_e);
pdfConII=connec_prob(sig(II),gridSize,R_i);
pdfConIE=connec_prob(sig(IE),gridSize,R_i,R_e);
pdfConEI=connec_prob(sig(EI),gridSize,R_e,R_i);

pdfConEE2=connec_prob(sig(EE2), gridSize, R_e, R_e2);

pdfConE2E2=connec_prob(sig(E2E2), gridSize, R_e2, R_e2);
pdfConE2I2=connec_prob(sig(E2I2), gridSize, R_e2, R_i2);
pdfConI2I2=connec_prob(sig(I2I2), gridSize, R_i2, R_i2);
pdfConI2E2=connec_prob(sig(I2E2), gridSize, R_i2, R_e2);

%% generate adjacency matrices
[A_IE]=sampleConn(pdfConIE,p(IE));
[A_EE]=sampleConn(pdfConEE,p(EE));
[A_II]=sampleConn(pdfConII,p(II));
[A_EI]=sampleConn(pdfConEI,p(EI));

[A_EE2]=sampleConn(pdfConEE2,p(EE2));
[A_E2E2]=sampleConn(pdfConE2E2,p(E2E2));
[A_E2I2]=sampleConn(pdfConE2I2,p(E2I2));
[A_I2I2]=sampleConn(pdfConI2I2,p(I2I2));
[A_I2E2]=sampleConn(pdfConI2E2,p(I2E2));




%% transform to synaptic weights
connMat_sep=cell(9,1);
connMat_sep{EE}=bsxfun(@rdivide, A_EE, sum(A_EE,2))*strength(EE);
connMat_sep{II}=bsxfun(@rdivide, A_II, sum(A_II,2))*strength(II);
connMat_sep{IE}=bsxfun(@rdivide, A_IE, sum(A_IE,2))*strength(IE);
connMat_sep{EI}=bsxfun(@rdivide, A_EI, sum(A_EI,2))*strength(EI);

connMat_sep{EE2}=bsxfun(@rdivide, A_EE2, sum(A_EE2,2))*strength(EE2);

connMat_sep{E2E2}=bsxfun(@rdivide, A_E2E2, sum(A_E2E2,2))*strength(E2E2);
connMat_sep{I2I2}=bsxfun(@rdivide, A_I2I2, sum(A_I2I2,2))*strength(I2I2);
connMat_sep{I2E2}=bsxfun(@rdivide, A_I2E2, sum(A_I2E2,2))*strength(I2E2);
connMat_sep{E2I2}=bsxfun(@rdivide, A_E2I2, sum(A_E2I2,2))*strength(E2I2);


connMat =zeros(totalN);
connMat(E_ind,E_ind)= connMat_sep{EE};
connMat(I_ind,E_ind)= connMat_sep{EI};
connMat(E_ind,I_ind)= connMat_sep{IE};
connMat(I_ind,I_ind)= connMat_sep{II};

connMat(E2_ind,E_ind)= connMat_sep{EE2};

connMat(E2_ind,E2_ind)= connMat_sep{E2E2};
connMat(I2_ind,E2_ind)= connMat_sep{E2I2};
connMat(E2_ind,I2_ind)= connMat_sep{I2E2};
connMat(I2_ind,I2_ind)= connMat_sep{I2I2};


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
DX=min(DX,gridSize(1)-DX); % to make it circular
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
cpdf=cumsum(pdfCon,2);% make comulative distr.
%nmber of incoming connections per neuron (actually number of samples, therefore maximum)
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

