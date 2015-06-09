function [output,spikes]=Izh_network(V_init,I,noise,cfg)
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
tLim=cfg.tLim;
dt=cfg.dt;

t=tLim(1):dt:tLim(2);

output.t=t;

numIt=numel(t);
numNeur=size(I,2)-1;


I_inp=zeros(numIt,numNeur);

% interpolate I and noise to desired resolution
% Note: perhaps change this to do this in steps to save memory space
I_orig=I;
I=interp1(I(:,1),I(:,2:end),t);

if nargin<3
  noise=sparse(numIt,numNeur);
else
  noise=interp1(noise(:,1),noise(:,2:end),t);
end

V=V_init;

%% parsing rest of cfg
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

try
  STDPflag=cfg.STDP;
catch
  STDPflag=true;
end


if STDPflag
  try
    tau_STDP=cfg.tau_STDP(:);
  catch    
    tau_STDP=[15; 15];
  end
  try
    A_STDP=cfg.A_STDP(:);
  catch
    A_STDP=[.01; .01];
  end
  
  E_syn=false(size(S));
  E_syn(:,1:sum(EI))=true;
  
  deltaS=zeros(numIt,2);
end

try
  verboseFlag=cfg.verbose;
catch
  verboseFlag=false;
end

try
  outpFname=cfg.output;
  outpFlag=true;
  
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
catch
  outpFlag=false;
end


if outpFlag
  try
    saveInterval=cfg.saveInterval;
  catch
    saveInterval=250;
  end  
end

try
  fullOutput=cfg.fullOutput;
catch
  fullOutput=0;
end

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

%% integration loop
for n=2:numIt
  
  I_E=G(1,EI)*S(:,EI).'.*(V_AMPA-V(1,:));
  I_I=G(1,~EI)*S(:,~EI).'.*(V_GABA-V(1,:));
  I_tot=I(n,:)+I_E+I_I;

  I_tot=I_tot+noise(n,:).*randn(1,numNeur);
  
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
    G(1,firSel)=1;
    
    if STDPflag
      if any(firSel)% update synapses using STDP
        % update STDP memory (Only E-E and E-I interactions)
        X(:,firSel)=bsxfun(@plus,X(:,firSel),A_STDP*0+1);
        
        if fullOutput
          X_list(:,:,n)=X(:,:);
        end
        
        dsyn=STDP(t,S,X(:,:),A_STDP,firSel,S>0 & E_syn);
        S=S+dsyn;
        deltaS(n,1)=sum(dsyn(:));
        deltaS(n,2)=sum(abs(dsyn(:)));
        % hard bound to be larger than 0
        S(S(:)<0)=0;
        % mulitplicative normalization
        sel=sum(S(EI,EI),2)>Smax;
        if any(sel)
          S(EIind(sel),EI)=bsxfun(@times,S(EIind(sel),EI),1./sum(S(EIind(sel),EI),2))*Smax;
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
  
  % save every 250 iterations (but skip the last one, if simulation is
  % almost done)
  if outpFlag && mod(n,saveInterval)==0 && (numIt-n)>saveInterval/2
    output.S=S;
    output.spiks=spiks;
    
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
    save(outpFname,'output')
    try
      save(outpFname,'spikes','-append')
    end
  end
  
  if verboseFlag
    msg=sprintf(['Iteration %d/%d\n'], [n numIt]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
  end
  
end


output.S=S;
output.spiks=spiks;

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
if outputFlag
  save(outpFname,'output')
  try
    save(outpFname,'spikes','-append')
  end
end
