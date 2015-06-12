function dsyn=STDP_izh2004(dt,firSel,lastSpike,A_stdp,tau_stdp,S_struct)
% dsyn=STDP_izh2004(firSel,spiks,S_struct)
% 
% Izhikevich EM, Gally J a., Edelman GM. Spike-timing dynamics of neuronal groups. Cereb Cortex. 2004
% 
% lastSpike should contain number of time samples since last spike of the neurons
% (will be converted using dt)
% 
% A_stdp = [A-, A+]
% tau_stdp=[tau-, tau+]


if nargin < 6
  S_struct=logical(size(spiks,1));
end

if iscolumn(S_struct) || isrow(S_struct)
  S_struct=reshape(S_struct,sqrt(numel(S_struct)),sqrt(numel(S_struct)));
end

if islogical(firSel)
  firIdx=find(firSel);
else
  firIdx=firSel;
end

dsyn=zeros(size(lastSpike,1));

for n=1:numel(firIdx)
  % Update as presynaptic neurons
  % presyn neuron j spikes, s'_ij decreases by A_stdp(1)*exp(t_i-t)/tau_stdp(1)
  postsyns=find(S_struct(n,:));
  delta_t=lastSpike(postsyns)/dt;
  dsyn(postsyns,firIdx(n))=-A_stdp(1)*exp(-delta_t./tau_stdp(1));
  
  % Update as postsynaptic neurons
  % postsyn neuron i spikes, s'_ij increases by A_stdp(2)*exp(t_j-t)/tau_stdp(2)
  presyns=find(S_struct(:,n));  
  delta_t=lastSpike(presyns)/dt;
  dsyn(firIdx(n),presyns)=A_stdp(2)*exp(-delta_t./tau_stdp(2));
  
end