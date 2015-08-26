function [ dS ] = TAH(t,S,X,A,pars,firSel,S_struct)
% [ dS ] = TAH(t,S,X,A,firSel,S_struct)
% 
% Dynamics for STDP (TAH; temporally asymmetric hebbian) , using STDP memory trace X
% from : 
% Gütig R, Aharonov R, Rotter S, Sompolinsky H. 
% Learning input correlations through nonlinear temporally asymmetric Hebbian plasticity. 
% J Neurosci. 2003;23: 3697–3714. 
% 
% Assumes S contains normalized synatpic couplings ([0 1])
% pars contains both exponent mu and asymmetry parameter alpha: [mu alpha];
% 
% Used for calling from RK4.

A_pos=A(1);
A_neg=A(2);

dS = zeros(size(S));

for firNeur=find(firSel)
  if nargin>5
    pre=S_struct(firNeur,:);
    post=S_struct(:,firNeur);
  else
    pre=true(1,size(S,2));
    post=false(size(S,1),1);
  end
  if any(pre)
    dS(firNeur,pre)=A_pos*X(1,pre).*(1-S(firNeur,pre)).^pars(1);
  end
  if any(post)
    dS(post,firNeur)=-A_neg*X(2,post).*S(post,firNeur).'.^pars(1)*pars(2);
  end
end
