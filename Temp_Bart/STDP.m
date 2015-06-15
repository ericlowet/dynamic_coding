function [ dS ] = STDP(t,S,X,A,firSel,S_struct)
% [ dS ] = STDP(t,S,X,A,firSel,S_struct)
% 
% Dynamics for STDP, using STDP memory trace X
% from : http://www.scholarpedia.org/article/STDP
% Used for calling from RK4.

A_pos=A(1);
A_neg=A(2);

dS = zeros(size(S));

% for pre=find(firSel)
%   for post=find(firSel)
%     if nargin<6 || S_struct(post,pre)
%       dS(post,pre)=firSel(post)*A_pos*X(1,pre)-firSel(pre)*A_neg*X(2,post);
%     end
%   end
% end

for firNeur=find(firSel)
  if nargin>5
    pre=S_struct(firNeur,:);
    post=S_struct(:,firNeur);
  else
    pre=true(1,size(S,2));
    post=false(size(S,1),1);
  end
  dS(firNeur,pre)=A_pos*X(1,pre);
  dS(post,firNeur)=-A_neg*X(2,post);
end
