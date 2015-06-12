function [spec, f]=bg_fftwelch(dat,fs,tapLen,demean,taper,nfft)
% [spec, f]=bg_fftwelch(dat,fs,tapLen,demean,taper,nfft)
% 
% dat     = data on which to calculate fft (timepoints x other dimensions)
% fs      = sampling rate of dat
% tapLen  = length of tapers/windows (in seconds); leaving tapLen empty or
%           setting it to 0 will result in using only 1 windows (i.e. 
%           performing a simple fft on the data)
% demean  = (optional; default=0) demean the windows, before calculating fft
% taper   = (optional; default=0) apply hamming taper to windows
% nfft    = (optional; default=0) same as nfft input of matlab's fft
%           function. Can be used to apply zero-padding to yield
%           interpolated spectra.


if isrow(dat)
  dat=dat(:);
end

if nargin<4
  demean=0;
end
if demean
  if iscolumn(dat)
    dat=dat-mean(dat);
  else
    dat=bsxfun(@minus,dat,mean(dat,1));
  end
end


datLen=size(dat,1);
if nargin<3 || isempty(tapLen) || tapLen==0 || (tapLen*fs)>size(dat,1)
  L=size(dat,1);
else
  L=ceil(tapLen*fs);
  dat=bg_reshape_overlap(dat,ceil(.5*L),L,1);
end


if nargin>4 && taper
  tap=hamming(L);
  dat=bsxfun(@times,dat,tap);
end

if nargin<6
  nfft=L;
elseif nfft<=1
  nfft=round(datLen*nfft);
end


df=1/nfft*fs;
if mod(nfft,2)
  f=[0:df:fs/2 -fs/2:df:0-df]';
else
  f=[0:df:fs/2-df -fs/2:df:0-df]';
end

spec=fft(dat,nfft)/(L);
% 
% f=fftshift(f);
% spec=fftshift(spec);