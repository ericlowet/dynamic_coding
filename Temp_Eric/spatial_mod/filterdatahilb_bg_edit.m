function [sig] = filterdatahilb_bg_edit(signals, LowF,HihF,fs)


[y] = ft_preproc_bandpassfilter(signals, fs, [LowF HihF ]);

sig= hilbert(y.').';
end














