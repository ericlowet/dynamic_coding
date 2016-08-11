function  [ft] = fft_dpss_sig(lfp)

Fs=1000;
L=length(lfp);
y=(lfp);   
wins = dpss(length(y),9);
for tapl=1:size(wins,2)
y=bsxfun(@times,y,wins(:,tapl)');
Y = 2.*(fft(y,L,2)./L);
f = Fs/2*linspace(0,1,L/2+1);
ft.fourierspctrm(tapl,:) = Y(:,1:L/2+1);
end
ft.fourierspctrm=mean(abs(ft.fourierspctrm));
ft.freq=f;