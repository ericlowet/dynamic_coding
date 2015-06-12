function [V,t,output,spikes]=Izh_network_STDP_izh2004(V_init,tLim,dt,I,noise,cfg)
% [V,t,output,spikes]=Izh_network(V_init,tLim,dt,I,noise,cfg)
% 
% simulation parameters:
% tLim  = [tStart tEnd]
% dt    = length of timestep
% 
% input to neurons:
% V_init= initial membrane potential for every neuron
% I   = input current to the neuron
% noise= standard deviation of gaussian white noise added to I.
% 
% cfg contains fields:
% a-d = parameters for the Izhikevich model
% S   = connection matrix containing maximum conductance values
% EI  = neuron label; determinse whether the neuron is excitatory or not 
%       (1= E; 0= I)
% 
% OPTIONAL fields:
% STDP    = flag to perform (or omit) spike timing-dependent plasticity
% (default=0; no STDP)
% verbose = flag to print progress (iteration number) to screen (default
% =1)
% 

%% Initializing variables
t=tLim(1):dt:tLim(2);
numIt=numel(t);
numNeur=size(I,2)-1;
V=nan(2,numNeur,numIt);
I_inp=zeros(numIt,numNeur);

% interpolate I and noise to desired resolution
I_orig=I;
I=interp1(I(:,1),I(:,2:end),t);

if nargin<5
  noise=sparse(numIt,numNeur);
else
  noise=interp1(noise(:,1),noise(:,2:end),t);
end

V(1,:,1)=V_init(1,:);
V(2,:,1)=V_init(2,:);

%% parsing cfg
try
  a=cfg.a;
  b=cfg.b;
  c=cfg.c;
  d=cfg.d;
catch
  error('cfg does not contain all paremeters for neurons (cfg.a-d)')
end

try
  S=cfg.S;  
catch
  warning('No synaptic connections given; cfg.S missing')
end


try
  EI=cfg.EI;
catch
  warning('Character of neurons undefined (cfg.EI); assuming all neurons are excitatory')
  EI=true(size(S,1));
end

output=[];
output.I=I_orig;
output.EI=EI;

try
  STDPflag=cfg.STDP;
catch
  STDPflag=true;  
end

if STDPflag
  try
    tau_STDP=cfg.tau_STDP;
  catch    
    tau_STDP=[20; 15];
  end
  try
    A_STDP=cfg.A_STDP;
  catch
    A_STDP=[.004; .004];
  end
end


try
  outpFname=cfg.output;
  outpFlag=true;
catch
  outpFlag=false;
end

try
  verboseFlag=cfg.verbose;
catch
  verboseFlag=false;
end

if outpFlag
  try
    saveInterval=cfg.saveInterval;
  catch
    saveInterval=250;
  end
end

%%


% conductances; AMPA and GABA -- NOTE this is FROM the neuron
tau_G=[10; 5];
G=nan(1,numNeur,numIt);
G(:,:,1)=0;



Smax=mean(sum(S(EI,EI),2));

spiks=sparse(numNeur,numIt);

V_AMPA=0;
V_GABA=-90;

if nargout>3
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

if STDPflag
  output.S_orig=S;  
  S_structSQ=logical(S);
  EISel=false(numNeur,numNeur);
  EISel(1:sum(EI),1:sum(EI))=true;
  
  lastSpike=ones(sum(EI),1)*inf;
  
  % keep track of dynamics of synapses that exist structurally and are E->E
  S_dyn=full([S(S_structSQ & EISel)'; 0*S(S_structSQ & EISel)']); % first temporal derivitive equal to 0
  A=1e-6;
end

%% integration loop
for n=2:numIt
  
  I_E=G(1,EI,n-1)*S(:,EI).'.*(V_AMPA-V(1,:,n-1));
  I_I=G(1,~EI,n-1)*S(:,~EI).'.*(V_GABA-V(1,:,n-1));
  I_tot=I(n,:)+I_E+I_I;

  I_tot=I_tot+noise(n,:).*randn(1,numNeur);
  I_inp(n,:)=I_tot;
  
  % membrane potential
  V(:,:,n)=RK4(t(n-1),V(:,:,n-1),dt,'Izh_neuron',a,b,I_tot);
  
  % synaptic conductance
  G(:,EI,n)=RK4(t(n-1),G(:,EI,n-1),dt,'exp_decay',tau_G(1));
  G(:,~EI,n)=RK4(t(n-1),G(:,~EI,n-1),dt,'exp_decay',tau_G(2));
  
  if STDPflag
    % Synaptic dynamics (Only E-E interactions)
    S_dyn(:,:,n)=RK4(t(n-1),S_dyn(:,:,n-1),dt,'SynDecay_izh2004',A);
    
    %manually clip synapses to [0 .05];
    synLims=[0 .05];
    S_dyn(1,S_dyn(1,:,n)<synLims(1),n)=synLims(1);
    S_dyn(1,S_dyn(1,:,n)>synLims(2),n)=synLims(2);
    
    S(S_structSQ & EISel)=S_dyn(1,:,n);
  end
    
  
  
  firSel=squeeze(V(1,:,n)>30);
  firSel_E=[];
  if any(firSel)
    
    % reset membrane potential
    if numel(c)>1
      V(1,firSel,n)=c(firSel);
      V(2,firSel,n)=V(2,firSel,n-1)+d(firSel);
    else
      V(1,firSel,n)=c;
      V(2,firSel,n)=V(2,firSel,n-1)+d;
    end
    
    % update synaptic channels to fully open
    G(1,firSel,n)=1;
    
    spiks(firSel,n)=true;
    % create spikes structure
    if nargout>3
      firSel=find(firSel);
      for firIdx=1:numel(firSel)
        spikes.timestamp{firSel(firIdx)}=[spikes.timestamp{firSel(firIdx)} t(n)];
      end
    end
    
    if STDPflag && n>2      
      firSel_E=firSel(firSel<=sum(EI));
      if ~isempty(firSel_E)% update synapses using STDP
        % update temporal derivative of synapses (Only E-E interactions)
        dsyn=STDP_izh2004(dt,firSel_E,lastSpike,A_STDP,tau_STDP,S_structSQ(1:sum(EI),1:sum(EI)));
        S_dyn(2,:,n)=S_dyn(2,:,n)+dsyn(S_structSQ(EISel))';
      end
    end
    
  end
  
  % save every 'saveInterval' iterations (but skip the last one, if simulation is
  % almost done)
  if outpFlag && mod(n,saveInterval)==0 && (numIt-n)>saveInterval/2
    output.G=G;
    output.S=S;
    output.spiks=spiks;
    if STDPflag
      output.S_dyn=S_dyn;
    end
    save(outpFname,'V','t','output')
    try 
      save(outpFname,'spikes','-append')
    end
  end
  
  if verboseFlag
    msg=sprintf(['Iteration %d/%d\n'], [n numIt]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  
  if STDPflag 
    lastSpike(firSel_E)=0;
    lastSpike=lastSpike+1;
  end
  
end

% output.I_tot=I_tot;
output.G=G;
output.S=S;
output.spiks=spiks;
if STDPflag
  output.S_dyn=S_dyn;
end
if outpFlag
  save(outpFname,'V','t','output')
  try
    save(outpFname,'spikes','-append')
  end
end

