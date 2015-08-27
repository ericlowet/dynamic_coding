function [output,spikes]=Izh_network_TAH(input)
% [V,t,output,spikes]=Izh_network_TAH(input)
% 
% 
% 
% input contains fields:
% a-d = parameters for the Izhikevich model
% S   = connection matrix containing maximum conductance values
% EI  = neuron label; determinse whether the neuron is excitatory or not 
%       (1= E; 0= I)
% 
% V_init= initial membrane potential for every neuron
% I   = input current to the neuron
% noise= standard deviation of gaussian white noise added to I.
% 
% tLim  = [tStart tEnd]
% dt    = length of timestep
% 
% OPTIONAL fields:
% STDP    = flag to perform (or omit) spike timing-dependent plasticity
% (default=0; no STDP)
% verbose = flag to print progress (iteration number) to screen (default
% =1)
% 
% Dynamics for STDP (TAH; temporally asymmetric hebbian) , using STDP memory trace X
% from : 
% Gütig R, Aharonov R, Rotter S, Sompolinsky H. 
% Learning input correlations through nonlinear temporally asymmetric Hebbian plasticity. 
% J Neurosci. 2003;23: 3697–3714. 

%% Initializing variables
tic
tLim=input.tLim;
dt=input.dt;

t=tLim(1):dt:tLim(2);


try
  I=input.I;
catch
  error('neurons have no input (input.I)')
end

try 
  noise=input.noise;
catch
  noise=zeros(2,numNeur);
  noise=[tLim(:) noise];
end

numIt=numel(t);
numNeur=size(I,2)-1;

I_inp=zeros(numIt,numNeur);

%% parsing rest of input
try
  a=input.a;
  b=input.b;
  c=input.c;
  d=input.d;
catch
  error('input does not contain all parameters for neurons (input.a-d)')
end

try
  V_init=input.V_init;
catch
  V_init=ones(2,numNeur);
  V_init=bsxfun(@times,V_init,[-70;0]); % -70 mV
  input.V_init=V_init;
end

V=V_init;

try
  S=input.S;
catch
  warning('No synaptic connections given; input.S missing')
end

try
  EI=input.EI;
catch
  warning('Character of neurons undefined (input.EI); assuming all neurons are excitatory')
  EI=true(size(S,1));
  input.EI=EI;
end

try
  STDPflag=input.STDP;
catch
  STDPflag=false;
end


if STDPflag
  try
    tau_STDP=input.tau_STDP(:);
  catch    
    tau_STDP=[20; 20];
    input.tau_STDP=tau_STDP;
  end
  try
    A_STDP=input.A_STDP(:);
  catch
    A_STDP=[.001; .001];
    input.A_STDP=A_STDP;
  end
  
  E_syn=false(size(S));
  E_syn(1:sum(EI),1:sum(EI))=true;
  
  deltaS=zeros(numIt,2);
end

try
  verboseFlag=input.verbose;
catch
  verboseFlag=true;
end

try
  outpFname=input.output;
  if ~isempty(outpFname)
    outpFlag=true;
    
    [dirN,fileN,ext]=fileparts(outpFname);
    if ~isempty(dirN) && ~exist(dirN,'file')
      mkdir(dirN)
    end
    if isempty(ext)
      % force .mat extension
      outpFname=fullfile(dirN,fileN,'.mat');
    end
    
    if verboseFlag && exist([outpFname],'file')
      warning([outpFname ' already exists...'])
      reply = input('Do you want to overwrite? Y/N [N]: ', 's');
      if isempty(reply)
        outpFlag = false;
      elseif strcmpi(reply,'n')
        outpFlag = false;
      elseif strcmpi(reply,'y')
        outpFlag = true;
      else
        warning('input not recognized, not overwriting')
        outpFlag = false;
      end
      
    end
    
    
    
  else
    outpFlag=false;
  end
catch
  outpFlag=false;
end


if outpFlag
  try
    saveInterval=input.saveInterval;
  catch
    saveInterval=250;
  end  
end

try
  fullOutput=input.fullOutput;
catch
  fullOutput=0;
end

try 
  TAHpars=input.TAHpars;
catch
  TAHpars=[.5, 1];
  input.TAHpars=TAHpars;
end

if numel(TAHpars)==1
  TAHpars=[TAHpars 1];
  input.TAHpars=TAHpars;
end

try
  maxSynVal=input.maxSynVal;
catch
  maxSynVal=max(max(S(EI,EI)));
  input.maxSynVal=maxSynVal;
end
  
try
  delay=input.delay; % in ms;
catch
  delay=0;
  input.delay=0;
end

if delay
  delayFlag=1;
else
  delayFlag=0;
end

%save input structure
output.input=input;


%%

% conductances; AMPA and GABA -- NOTE this is FROM the neuron
tau_G=[10; 5];
G=zeros(1,numNeur);

% STDP memory
if STDPflag
  X=zeros(2,numNeur);
end

Smax=mean(sum(S(EI,EI),2));

spiks=sparse(numNeur,numIt);

V_AMPA=50;
V_GABA=-90;

if nargout>1
  spikes=[];
  spikes.label=cell(1,numNeur);
  spikes.label(EI)={'E'};
  spikes.label(~EI)={'I'};
  spikes.timestamp=cell(1,numNeur);
end

EIind=find(EI);

if verboseFlag
  reverseStr=[];
end

I_orig=I;
noise_orig=noise;
interpIval=1e3;

output.S_orig=S;
output.I=I_orig;
output.EI=EI;

if fullOutput
  V_list=nan(2,numNeur,numIt);  
  G_list=nan(1,numNeur,numIt);
  G_list(:,:,1)=0;
  if STDPflag
    X_list=nan(2,numNeur,numIt);
    X_list(:,:,1)=0;
  end
end



% interpolate I and noise to desired resolution
% Note: this is done in steps of interpIval to save memory space
I=interp1(I_orig(:,1),I_orig(:,2:end),t([1:interpIval]));
noise=interp1(noise_orig(:,1),noise_orig(:,2:end),t([1:interpIval]));
I(isnan(I))=0;
noise(isnan(noise))=0;
interpCnt=1;

if delayFlag
  firSelHist=false(delay/dt+1,numNeur);
end

%% integration loop
for n=2:numIt
  
  interpCnt=interpCnt+1;
  if interpCnt==(interpIval+1)
    % interpolate I and noise to desired resolution
    % Note: this is done in steps of interpIval to save memory space
    interpCnt=1;
    tVec=[1:interpIval]+n-1;
    tVec=tVec(tVec<=numel(t));
    I=interp1(I_orig(:,1),I_orig(:,2:end),t(tVec));
    noise=interp1(noise_orig(:,1),noise_orig(:,2:end),t(tVec));
    I(isnan(I))=0;
    noise(isnan(noise))=0;
  end
  
  I_E=G(1,EI)*S(:,EI).'.*(V_AMPA-V(1,:));
  I_I=G(1,~EI)*S(:,~EI).'.*(V_GABA-V(1,:));
  I_tot=I(interpCnt,:)+I_E+I_I;

  I_tot=I_tot+noise(interpCnt,:).*randn(1,numNeur);
  
  % membrane potential
  V(:,:)=RK4(t(n-1),V(:,:),dt,'Izh_neuron',a,b,I_tot);
    
  % synaptic conductance
  G(:,EI)=RK4(t(n-1),G(:,EI),dt,'exp_decay',tau_G(1));
  G(:,~EI)=RK4(t(n-1),G(:,~EI),dt,'exp_decay',tau_G(2));

  if STDPflag
    % STDP memory
    X(:,:)=RK4(t(n-1),X(:,:),dt,'exp_decay',tau_STDP);
    if fullOutput
      X_list(:,:,n)=X(:,:);
    end
  end
    
  
  
  firSel=squeeze(V(1,:)>30);
  
  if delayFlag
    firSelHist(mod(n,delay/dt+1)+1,:)=firSel;
  end
  
  if any(firSel)
    
    % reset membrane potential
    if numel(c)>1
      V(1,firSel)=c(firSel);
      V(2,firSel)=V(2,firSel)+d(firSel);
    else
      V(1,firSel)=c;
      V(2,firSel)=V(2,firSel)+d;
    end
    
    % update synaptic channels to fully open
    if delayFlag
      G(1,firSelHist(mod(n+1,delay/dt+1)+1,:))=1;
    else
      G(1,firSel)=1;
    end
    
    if STDPflag
      if any(firSel)% update synapses using STDP
        dsyn=TAH(t,S/maxSynVal,X(:,:),A_STDP,TAHpars,firSel,S>0 & E_syn);
        dsyn=dsyn*maxSynVal; % synapses are bounded/normalized by maxSynVal
        dsyn(dsyn<0)=max(dsyn(dsyn<0),-S(dsyn<0)); % clip, because negative values will yield synapses with imaginary components
        S=S+dsyn;
        deltaS(n,1)=sum(dsyn(:));
        deltaS(n,2)=sum(abs(dsyn(:)));
        
        % update STDP memory (Only E-E and E-I interactions)
        % but only after doing STDP, cells firing at exactly the same time
        % have no influence        
%         X(:,firSel)=bsxfun(@plus,X(:,firSel),A_STDP*0+1);
        X(:,firSel)=1;
        
        if fullOutput
          X_list(:,:,n)=X(:,:);
        end
        
      end
    end

    spiks(firSel,n)=true;
    % create spikes structure
    if nargout>1
      firSel=find(firSel);
      for firIdx=1:numel(firSel)
        spikes.timestamp{firSel(firIdx)}=[spikes.timestamp{firSel(firIdx)} t(n)];
      end
    end
    
  end
  
  % save every couple of iterations (but skip the last one, if simulation is
  % almost done)
  if outpFlag && mod(n,saveInterval)==0 && (numIt-n)>saveInterval/2
    output.S=S;
    output.spiks=spiks(:,1:n);
    output.t=t(1:n);
    output.simulationTime=toc;
    
    if STDPflag
      output.deltaS=deltaS(1:n,:);
      if fullOutput
        output.X=X_list;
      end
    end    
    
    if fullOutput
      G_list(:,EI,n)=G(:,EI);
      G_list(:,~EI,n)=G(:,~EI);
      V_list(:,:,n)=V;
      I_inp(n,:)=I_tot;
      output.G=G_list(:,:,1:n);
      output.I_tot=I_inp(1:n,:);
      output.V=V_list(:,:,1:n);
    end
    save(outpFname,'output','-v7.3')
    try
      save(outpFname,'spikes','-append','-v7.3')
    end
  end
  
  if verboseFlag
    msg=sprintf(['Iteration %d/%d\n'], [n numIt]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  
end

output.t=t;
output.S=S;
output.spiks=spiks;
output.simulationTime=toc;

if STDPflag
  output.deltaS=deltaS;
  if fullOutput
    output.X=X_list;
  end
end

if fullOutput
  G_list(:,EI,n)=G(:,EI);
  G_list(:,~EI,n)=G(:,~EI);
  V_list(:,:,n)=V;
  I_inp(n,:)=I_tot;
  output.G=G_list;
  output.I_tot=I_inp;
  output.V=V_list;
end
if outpFlag
  save(outpFname,'output','-v7.3')
  try
    save(outpFname,'spikes','-append','-v7.3')
  end
end
