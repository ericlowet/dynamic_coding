a=[0.02*ones(Ne,1);        0.1*ones(Ni,1)]'; % a =tiemscale of recover varibale u
b=[0.2*ones(Ne,1);         0.2*ones(Ni,1)]'; % b= sensitivity of  u to subthreshold oscillations
c=[-65*ones(Ne,1);        -65*ones(Ni,1)]; % c= membrane voltage after spike (reset)
d=[8*ones(Ne,1);           2*ones(Ni,1)]'; % d= spike reset of recover varibale u
v=-65*ones(Ne+Ni,1);    % Initial values of v = voltage
u=b'.*v;                % Initial values of u= membrane recovery variable

tau_G=[2; 8];
V_AMPA=50;V_GABA=-90;
G=zeros(1,Ntot,'single'); spiks=zeros(Ntot,simulation_time/dt,'single');
V(1,:)= ones(1,Ntot).*-70+randn(1,Ntot).*10;
V(2,:)= zeros(1,Ntot).*0+randn(1,Ntot).*0.1;
E1_ind=1:Ne; Etot=1:Ne;
I1_ind=Ne+1:Ntot;Itot=I1_ind;
firings=[];             % spike timings