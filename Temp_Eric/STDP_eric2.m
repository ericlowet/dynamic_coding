function [ dS ] = STDP_eric2(t,S,X,A,firSel,S_struct)
% [ dS ] = STDP(t,S,X,A,firSel,S_struct)
%
% Dynamics for STDP, using STDP memory trace X
% from : http://www.scholarpedia.org/article/STDP
% Used for calling from RK4.

A_pos= A(1);
A_neg=-A(2);
dS = zeros(size(S));
for spike_n=find(firSel)
    % look on the memory variable X of all neurons conected to spiking
    % neuron
    dS(spike_n, S_struct(spike_n,:)) = X(1,S_struct(spike_n,:)).*A_neg; % post-synaptic spike
    dS(S_struct(:,spike_n),spike_n)  = X(1,S_struct(:,spike_n)).*A_pos; % pre-synaptic spike
end
