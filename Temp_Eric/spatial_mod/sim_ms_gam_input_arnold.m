
%% spatial defined input pattern 
I2=sin(1.*linspace(-pi,pi,numNeurLin)+pi/2-0.0);
I2=repmat(I2,numNeurLin,1)';
I2=I2+I2';
I2(I2<0.4)=0;
I2(I2>0)=1;
%figure,imagesc(I2)

I=sin(2.*linspace(-pi,pi,numNeurLin)+pi/2-0.0);
I=repmat(I,numNeurLin,1);
I=randn(numNeurLin);%
%I=I+I';

sigI=4*dR;
Rdum=min(R_e):dR:max(R_e);
Rdum=Rdum-mean(Rdum);
filtKern=exp(-Rdum.^2/(2*sigI^2));
filtKern=bsxfun(@times,filtKern,filtKern');
filtKern=filtKern/sum(filtKern(:));
I=ifft2(fft2(I).*fft2(filtKern));
I=I-min(I(:));
I=((I/max(I(:))-0)*1).*mod_E;
%I=I.*I2;
%I(I<prctile(I(:),80))=0;
figure,imagesc(I)

I_2=I(1:2:end,1:2:end).*2/3;

%pause(1)
%%%%%%%%%%%%%%%%%%%
mean_E= [ E_inp1*ones(Ne1,1)'+I(:)' E_inp2*ones(Ne2,1)' ];
var_E1=4;var_E2= 4;
% to inhibitory neurons
mean_I= [ I_inp1*ones(Ni1,1)'+I_2(:)'  I_inp2*ones(Ni2,1)' ];
var_I= 4;
%%%%%%%%%%%

%%%%%%%%%
cfg=[];
cfg.simulation_time=simulation_time;
cfg.ms_interval=420;cfg.ms_interval_var=0;
cfg.sampling_frequency=2000;
cfg.modulation_strength=1;
cfg.positive_mod=0.8;cfg.negative_mod=0.4;
[ ms_signal,ms_times]  = create_ms_sig(cfg);
ms_signal=ms_signal(:).';
% figure,plot(ms_signal)%hold on,plot(ms_times-320,1,'o')
ms_dips=ms_times(2:end)-52;
ms_dips=ms_dips+round((randn(size(ms_dips)).*0));
ms_sig = single([ repmat(ms_signal,Ne1,1);repmat(ms_signal,Ne2,1).*0; repmat(ms_signal,Ni1,1).*1]);
ms_sig=ms_sig-min(ms_sig(:));