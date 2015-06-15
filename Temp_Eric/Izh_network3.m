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

%save input structure
output.cfg=cfg;

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
% parameter izhi neurons
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
%%%%%%%%%%%%%%%%%%
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
    A_STDP=[.0001; .0001];
  end
  
  E_syn=false(size(S));
  E_syn(1:sum(EI),1:sum(EI))=true;
  
  deltaS=zeros(numIt,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  outpFname=cfg.output;
  outpFlag=true;
  
%   if verboseFlag && exist([outpFname],'file')
%     warning([outpFname ' already exists...'])
%     reply = input('Do you want to overwrite? Y/N [N]: ', 's');
%     if isempty(reply)
%       outpFlag = false;
%     elseif strcmpi(reply,'n')
%       outpFlag = false;
%     elseif strcmpi(reply,'y')
%       outpFlag = true;
%     else
%       warning('input not recognized, not overwriting')
%       outpFlag = false;
%     end
%   end
catch
  outpFlag=false;
end
try
  fullOutput=cfg.fullOutput;
catch
  fullOutput=0;
end
try
  verboseFlag=cfg.verbose;
catch
  verboseFlag=true;
end

if outpFlag
  try
    saveInterval=cfg.saveInterval;
  catch
    saveInterval=250;
  end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


% G  =synaptic variable
% V  voltage and u
% 
%clear  firSelhist
firSelhist{1} =logical(zeros(length(EI),1,'int32'));
firSelhist{2} =firSelhist{1};
firSelhist{3} =firSelhist{1};
%% integration loop
for n=3:numIt
  n;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  I_E=G(1,EI)*S(:,EI).'.*(V_AMPA-V(1,:));
  I_I=G(1,~EI)*S(:,~EI).'.*(V_GABA-V(1,:));
  I_tot=I(n,:)+I_E+I_I;
  I_tot=I_tot+noise(n,:).*randn(1,numNeur);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % membrane potential
  V(:,:)=RK4(t(n-1),V(:,:),dt,'Izh_neuron',a,b,I_tot);
  % synaptic conductance
  G(:,EI)=RK4(t(n-1),G(:,EI),dt,'exp_decay',tau_G(1));
  G(:,~EI)=RK4(t(n-1),G(:,~EI),dt,'exp_decay',tau_G(2));

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  firSel=squeeze(V(1,:)>30);
  firSelhist{n-n+1} =firSel; 
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
    G(1,firSelhist{3})=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if STDPflag
    % STDP memory
    X(:,:)=RK4(t(n-1),X(:,:),dt,'exp_decay',tau_STDP);
    if fullOutput
      X_list(:,:,n)=X(:,:);
    end
  end
    
    if STDPflag
      if any(firSel)% update synapses using STDP
        % update STDP memory (Only E-E and E-I interactions)
        % Update X, but limit maximal level
       % X(:,firSel)=bsxfun(@plus,X(:,firSel),A_STDP*0+1);
       X(:,firSel)=1;
        
        
        if fullOutput
          X_list(:,:,n)=X(:,:);
        end
       
        
        dsyn=STDP_eric2(t,S,X(:,:),A_STDP,firSel,S>0 & E_syn);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    output.spiks=spiks(:,1:n);
    output.t=t(1:n);
    output.S_orig=cfg.S_orig;
    if STDPflag
      output.deltaS=deltaS(1:n,:);
      if fullOutput
        output.X=X_list;
      end
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
  
   firSelhist{3}= firSelhist{2} ;
   firSelhist{2}= firSelhist{1} ;

  
  
end   % loop 


output.t=t;
output.S=S;
output.S_orig=cfg.S_orig;
output.spiks=spiks;

if STDPflag
  output.deltaS=deltaS;
  if fullOutput
    output.X=X_list;
  end
end

%%%%% save in between data %%%%%%%%%
if outpFlag
  save(outpFname,'output')
  try
    save(outpFname,'spikes','-append')
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%