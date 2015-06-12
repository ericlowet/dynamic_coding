function [ dS ] = SynDecay_izh2004(t,S,a,dyn_label)
% [ dS ] = SynDecay_izh2004(t,S,a,S_struct)
% 
% dyn_label can be used to skip synapses that should stay stationary (only
% synapses mentioned in dyn_label are updated)
% 
% Dynamics for decaying synapses (2nd order differential eqn)
% S"=-(S'-a)*1e-4
% from
% Izhikevich EM, Gally J a., Edelman GM. Spike-timing dynamics of neuronal groups. Cereb Cortex. 2004
% Used for calling from RK4.



dS = zeros(size(S));

if nargin>3
  sel=dyn_label;
else
  sel=1:size(S,2);
end

dS(1,sel)=S(2,sel);
dS(2,sel)=-(S(2,sel)-a)*1e-4;
