%function sim_ms_gam_RF_batch29Juneb(iteration_sim)
clear all
addpath('/home/coherence/erilow/Network_modelling/')
addpath('/home/coherence/erilow/Network_modelling/sim_ms_RF/')
iteration_sim=1;


% vars1= [ 2 ];
%col =colormap(jet(length(vars1)));
f1=figure('Color','w')
%f2=figure('Color','w')
%f3=figure('Color','w')
nn=0;
%for VARI1=[vars1]
    nn=nn+1;
simulation_time=10020 ; %in ms
stim_size=0.5;stim_size2=2;
%% Overall input strenth
E_inp1=7;  E_inp2=1;I_inp1=E_inp1*2/3+1;I_inp2=2;
 mod_E=7.5;  %I_inp1A=E_inp1A*2/3;
clear spiks v

%% Lattice size 

numNeurLin=40; 

dR=50e-1;dt=0.5;  %time step

Ne1= numNeurLin^2;Ni1=Ne1/4;
Ne2=2^2;Ni2=Ne2/4;
Ne=Ne1+Ne2;  Ni=Ni1+Ni2;  
Ntot=(Ne+Ni);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_ms_gam_connectivity_arnold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_ms_gam_input_arnold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_ms_gam_init_AN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 sim_ms_gam_gen_Snoise
n=0;
for t=1:dt:simulation_time          % simulation
    n=n+1
    %    if mod(n,10)==0
     sim_ms_gam_gen_Snoise
  %  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_E=G(E1_ind)*S(:,E1_ind).'.*(V_AMPA-V(1,:))+[NoiseE+(mean_E.*ms_sig(1:Ne,n)') var_I.*randn(1,Ni)+(mean_I.*ms_sig(Ne:end,n)') ] ;
    I_I=G(Itot)*S(:,Itot).'.*(V_GABA-V(1,:));
    I_tot=I_E+I_I;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V(:,:)=RK4(t,V,dt,'Izh_neuron',a,b,I_tot);
    G(:,E1_ind)=RK4(t,G(:,E1_ind),dt,'exp_decay',tau_G(1));
    G(:,I1_ind)=RK4(t,G(:,I1_ind),dt,'exp_decay',tau_G(2));
    firSel=squeeze(V(1,:)>30);
    if any(firSel)
        if numel(c)>1
            V(1,firSel)=c(firSel);V(2,firSel)=V(2,firSel)+d(firSel);
        else
            V(1,firSel)=c; V(2,firSel)=V(2,firSel)+d;
        end
            G(1,firSel)=1;
        spiks(firSel,n)=1;
    end
   v1= V(1,:);v1(v1>-30)=-30;
    volt(:,n)= I_E.*1.2+I_I;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sp_raster_plot(int8(spiks(:,300:1:1400)),Ne1,Ne2,Ni1,Ni2)

%figure,imagesc(spiks)

%% Rates
St_rate=reshape(sum(spiks(1:Ne1,:),2),numNeurLin,numNeurLin)./(simulation_time./1000)   ;

figure('COlor','w','Position', [300 300 220 180])
h=subplot(1,1,1,'Fontsize',17)
imagesc(St_rate)
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
%colormap('jet')

All_rate=sum(spiks(1:Ntot,:),2)./(simulation_time./1000)   ;
%figure,plot(All_rate)
clear rates
rates.E1(1)=mean(All_rate(1:Ne1));rates.E1(2)=std(All_rate(1:Ne1));
rates.E2(1)=mean(All_rate(Ne1+1:Ne));rates.E2(2)=std(All_rate(Ne1+1:Ne));
rates.I1(1)=mean(All_rate(Ne+1:Ne+Ni1));rates.I1(2)=std(All_rate(Ne+1:Ne+Ni1));
rates.I2(1)=mean(All_rate(Ne+Ni1+1:Ntot));rates.I2(2)=std(All_rate(Ne+Ni1+1:Ntot));
spk_dens_E=fastsmooth(sum(spiks(1:Ne1,:)),60,3,1);


spk_dens_E2=fastsmooth(sum(spiks(Ne1+1:Ne,:)),60,3,1);

%% Signals
[ signals , selm] = make_sig_fast(volt(1:Ne1,:),R_e,numNeurLin);
% %signalV2=zscore(fastsmooth(sum(spiks(Ne1+1:Ne1+Ne2,1:2:end)),3,3,1));
%signals=volt(1:Ne1,1:2:end);
timl=size(signals,2);
[ ms_sp] = mstrigspk(ms_dips./1,spiks(1:Ne1,1:1:end));
t1=fastsmooth(mean(mean(ms_sp(selm,:,:),1),3),1,3,1);
allt1(:,nn)=t1;
t1=fastsmooth(t1,4,3,1);

[sigG] = filterdatahilb(signals, 22,65,1000);
xx=exp(i.*angle(sigG(:,1:1:end)))';vv=(xx' *xx);
x=round(numNeurLin/2);y=round(numNeurLin/2);
 seed= (x-1)*sqrt(numNeurLin^2)+y;
%%%%%%%%%%%%%%%%%%%%%%
%%% Input diff%%%%
clear allIdiff
for ind=1:Ne1
    for ind2=1:Ne1
        
        allIdiff(ind,ind2)= mean_E(ind)-mean_E(ind2);
    end
end

%%%%%%%%%%%%%%%%%
clear consum
for ind=1:Ne1
    for ind2=1:Ne1  
       consum(ind,ind2)= ([S3(ind,ind2)+S3(ind2,ind)])./2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%



[ ms_sigG] = mstrigsig(ms_dips./2,sigG);
[elnum,mstim, ms_n]=size(ms_sigG);
% %figure,plot( mean( mean(abs(ms_sigG(:,:,:)),1),3))
clear allarnos allarnos2 allarnosP allarnosP2
tt=randperm(1600);nn=0;
for seedC=[ tt(1:60)   ]   ;
    nn=nn+1
winti=5;nr=0;clear mscohs  msphs
for timsh=1:15:mstim-winti
    nr=nr+1;
    mGG=reshape(ms_sigG(:,timsh:timsh+winti,:),Ne1,(winti+1)*ms_n,1);
    mSP=reshape(ms_sp(:,timsh:timsh+winti,:),Ne1,(winti+1)*ms_n,1);
    
    xx=exp(1i.*angle(mGG))';vv=(xx' *xx);
    mscohs(:,:,nr)=single(reshape((abs(vv(seedC,1:Ne1))), numNeurLin, numNeurLin)./size(mGG,2));
    msphs(:,:,nr)=reshape((angle(vv(seedC,1:Ne1))), numNeurLin, numNeurLin);
   ms_timsel(nr)=([timsh+timsh+winti])./2;
%    XX=(exp(1i.*angle(mGG)));kl=corrcoef(XX');
% XX=double(angle(mGG));
% clear pp
% for ind=seedC
%     for ind2=1:1600     
%      pp(ind2)=  circ_corrcc(XX(ind,:),XX(ind2,:));        
%     end
% end
%     % mscohs2(:,:,nr)=  reshape(abs(kl(780,:)).*2,40,40);
%        mscohs2(:,:,nr)=  reshape(pp,40,40);
end
  selind= find(ms_timsel>=250 & ms_timsel <=350);
InpD=allIdiff(:,seedC);
ConS=consum(seedC,:);
Syn=(nanmean(mscohs(:,:,selind),3));
Phn=(circ_mean(msphs(:,:,selind),[],3));
 [ arno_map]=map_arnold(InpD(:), ConS(:), Syn(:));
 allarnos(:,:,nn)=arno_map;
  [ arno_map]=map_arnoldP(InpD(:), ConS(:), Phn(:));
 allarnosP(:,:,nn)=arno_map;
   selind= find(ms_timsel>=50 & ms_timsel <=90);
Syn=(nanmean(mscohs(:,:,selind),3));Phn=(nanmean(msphs(:,:,selind),3));
 [ arno_map]=map_arnold(InpD(:), ConS(:), Syn(:));
 allarnos2(:,:,nn)=arno_map;
   [ arno_map]=map_arnoldP(InpD(:), ConS(:), Phn(:));
 allarnosP2(:,:,nn)=arno_map; 
end
    [ arno_map, YY,XX]=map_arnold(InpD(:), ConS(:), Phn(:));
figure('COlor','w','Position', [300 300 240 180])
h=subplot(1,1,1,'Fontsize',17)
 imagesc(XX,YY,nanmean(allarnos,3))
axis xy;colormap('hot');set(gca,'Clim', [ 0.25 1])
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
figure('COlor','w','Position', [300 300 240 180])%%%%%%%%%%%%%%
h=subplot(1,1,1,'Fontsize',17)
 imagesc(XX,YY,nanmean(allarnos2,3))
axis xy;colormap('hot');
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
set(gca,'Clim', [ 0.25 1])
%%%%%%%%%%%%%%%%%%%%%%
figure('COlor','w','Position', [300 300 240 180])
h=subplot(1,1,1,'Fontsize',17)
 imagesc(XX,YY,circ_mean(allarnosP,[],3))
axis xy;colormap('hsv');set(gca,'Clim', [ -1 1])
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
figure('COlor','w','Position', [300 300 240 180])%%%%%%%%%%%%%
h=subplot(1,1,1,'Fontsize',17)
 imagesc(XX,YY,circ_mean(allarnosP2,[],3))
axis xy;colormap('hsv');%set(gca,'Clim', [ -1 1])
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[ ms_sig] = mstrigsigL(ms_dips./2,signals);
[elnum,mstim, ms_n]=size(ms_sigG);
figure('COlor','w','Position', [300 300 240 150])
h=subplot(1,1,1,'Fontsize',17)
plot((1:length(ms_sig(1,:,1)))./1000 -0.45, fastsmooth(ms_sig(780,:,7),1,3,1),'k')
set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
axis tight
xlim([ -0.08 0.46])



if 0
    clear data
    el=780;el2=80;
    tt=length(ms_sig(el,:,1));
    for tr=1:size(ms_sig,3);
        data.trial{tr}(1,:) =zscore(ms_sig(el,:,tr))+randn(1,tt).*1;data.trial{tr}(2,:) =zscore(ms_sig(el2,:,tr))+randn(1,tt).*1;
        data.time{tr} =(1:length(ms_sig(el,:,tr)))./1000;
        data.label(1)={'Channel1'};data.label(2)={'Channel2'};
    end
    cfg = [];
    cfg.method ='wavelet';%'mtmconvol' for STFT;
    cfg.output ='fourier';cfg.keeptrials ='yes'
    cfg.channel= 'all';cfg.foi= 8:1:65;
    cfg.toi= data.time{1}(1:1:end);
    cfg.width =6; % important parameter
    gh=pwd; cd('/home/common/matlab/fieldtrip/')
    waveTFR= ft_freqanalysis(cfg, data);
    cd(gh)
    k=figure('COlor','w','Position', [300 300 240 180])
    h=subplot(1,1,1,'Fontsize',17)
    imagesc(  waveTFR.time-0.45,  waveTFR.freq,squeeze(mean(abs(waveTFR.fourierspctrm(:,1,:,:))))    )
    axis xyxlim([ -0.08 0.4])
    set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
   
    pd= squeeze(circ_dist(angle(waveTFR.fourierspctrm(:,1,:,:)), angle(waveTFR.fourierspctrm(:,2,:,:)))  );
    p1= exp(1i.*pd);PLVTFR=squeeze(abs(nanmean(p1,1)));
    figure('COlor','w','Position', [300 300 240 180])
    h=subplot(1,1,1,'Fontsize',17)
    imagesc(waveTFR.time-0.45,  waveTFR.freq,PLVTFR)
    axis xy xlim([ -0.08 0.4])
    set(h,'FontName','Arial','FontSize',12,'FontWeight','bold');
    colormap('hot');set(gca,'Clim',[ 0 1])
end