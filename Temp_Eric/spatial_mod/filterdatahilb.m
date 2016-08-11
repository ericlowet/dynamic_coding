function [sig] = filterdatahilb(signals, LowF,HihF,fs)

[n,m]=size(signals);
sig=signals;
for el=1:n
[y] = ft_preproc_bandpassfilter(signals(el,:), fs, [LowF HihF ]);
sig(el,:)= hilbert(y);
end














