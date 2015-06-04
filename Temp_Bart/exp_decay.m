function [ dX ] = exp_decay(t,X,tau)
% [ dX ] = exp_decay(t,X,tau)
% 
% Dynamics for exponentially decaying variable (e.g. synaptic conductance)
% Used for calling from RK4.


% from 2008 paper
dX=-bsxfun(@rdivide,X,tau);