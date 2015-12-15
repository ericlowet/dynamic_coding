function [ out, ms_t]  = create_ms_sig(cfg)
% [ out, ms_t]  = create_ms_sig(cfg)
% Example usage:
% cfg=[];
% cfg.sampling_frequency=1000;
% cfg.simulation_time=5000 %milliseconds
% cfg.ms_interval=330
% cfg.ms_interval_var=30;
% cfg.modulation_strength=3;
% cfg.modulation_var=1e-2;
% cfg.positive_mod=1;
% cfg.negative_mod=0.25;
% [ ms_signal,ms_times]  = create_ms_sig(cfg);
%
% out=  ms-modulated input
% ms_t = time points of ms


if isfield(cfg,'simulation_time')
  t_max=cfg.simulation_time;     % simulation_time
else
  t_max =1000;
end

if isfield(cfg,'modulation_strength')
  max_mod =cfg.modulation_strength;     % ms_modulation strength
else
  max_mod =5;
end

if isfield(cfg,'ms_interval')
  interval =cfg.ms_interval;     % default ms interval
else
  interval=300;
end

if isfield(cfg,'ms_interval_var')
  msvar =cfg.ms_interval_var;     %  ms interval variability
else
  msvar=0;
end

if isfield(cfg,'modulation_var')
  modvar =cfg.modulation_var;     %  ms interval variability
else
  modvar=0;
end


if isfield(cfg,'sampling_frequency')
  FS =cfg.sampling_frequency;     % sampling_frequency in Hz
else
  FS=1000;
end

if isfield(cfg,'positive_mod')
  PM =cfg.positive_mod;     %  positive modulation peak in max_mod %  (default 1)
else
  PM=1;
end

if isfield(cfg,'negative_mod')
  NM =cfg.negative_mod;     % negative modulation peak in max_mod % (default 0.3)
else
  NM=0.3;
end

if isfield(cfg,'base_rate')
  base_rate=cfg.base_rate;
else
  base_rate=0;
end
  
%%%%%
interval=interval.*(FS/1000);
t_max=t_max.*FS/1e3;
msvar=msvar*FS/1e3;

baserate2=0; % not used

%% double exponentials %%%%%%
%%%% positive modulation part %%%%%
dt=1e3/FS;
timed=-450:dt:450;
C1=[ 0  0  +0 -60  1 ];
C2=[ 0  0  +0 -20  1 ];
[prob]= double_exp(timed,C1,C2,PM);
%%%% negative modulation part %%%%%
C1=[ 0  0  +0 -22  1 ];
C2=[ 0  0  +0 -15  1 ];
[prob2]= double_exp_rev(timed,C1,C2,-NM);
prob=prob+prob2;

z=(max_mod-baserate2).*(prob);

%% Create time series with ms times %%%%%%%%%

maxNumMs=ceil(t_max/interval*1.5);
ms=ones(maxNumMs,1)*interval+randn(maxNumMs,1)*msvar^.5;
ms_t=ceil(cumsum(ms));
ms_t=ms_t(ms_t<t_max);

out=zeros(t_max,1);
% add ms with variance in their amplitude. 
% (standard deviation divided by max_mod to get back to original units)
out(ms_t)=1+randn(numel(ms_t),1)*(modvar^.5/max_mod);


%%%%% Convolution %%%
out=conv(out,z,'same');
out=out+base_rate ;
%%%%%%%%%%%%%%%%%%%

