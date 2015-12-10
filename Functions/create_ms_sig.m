function [ out, ms_t]  = create_ms_sig(cfg) %t_max,max_mod,interval,msvar, FS)


% cfg=[];
% cfg.sampling_frequency=1000;
% cfg.simulation_time=5000
% cfg.ms_interval=330
% cfg.ms_interval_var=30;
% cfg.modulation_strength=3;
% cfg.positive_mod=1;
% cfg.negative_mod=0.25;
% [ ms_signal,ms_times]  = create_ms_sig(cfg);

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

clear prob u z
%%%%%
interval=interval.*(FS/1000);

baserate=0;
baserate2=0;

%% double exponentials %%%%%%
%%%% positive modulation part %%%%%
timed=-180:(1000/FS):450;
C1=[ 0  0  +0 -60  1 ];
C2=[ 0  0  +0 -20  1 ];
[prob]= double_exp(timed,C1,C2,PM);
%%%% negative modulation part %%%%%
timed=-180:(1000/FS):450;
C1=[ 0  0  +0 -22  1 ];
C2=[ 0  0  +0 -15  1 ];
[prob2]= double_exp_rev(timed,C1,C2,-NM);
prob=prob+prob2;


z=(max_mod-baserate2).*(prob)+baserate ;

%%%%%%%%%%%%%%%%
nt=1;
u=zeros(1,t_max);
u(cumsum(interval))=1;
u=u(1:t_max);
ms_t= cumsum(interval);
%%%%%%%%%%%%%%%%
%% Create time series with ms times %%%%%%%%%
if 1
    u=zeros(1,t_max);
    nx=1;
    ms_t=1;
    for t=1:length(u)-50
        if mod(t,interval(1)) ==0 & ms_t(nx)+(170./(1000./FS)) < t
            nx=nx+1;
            randshift= ceil(rand(1).*msvar);
            u(t+  randshift)=1;
            ms_t(nx)=t+randshift;
            nt=nt+1;
        end
    end
    % save('ms_t','ms_t')
end




%%%%% Convolution %%%
out=conv(u,z,'same');

%%%%%%%%%%%%%%%%%%%
 
