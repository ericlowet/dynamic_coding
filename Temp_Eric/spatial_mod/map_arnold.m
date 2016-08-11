function  [ arno_map, bin_cop, bin_freq]=map_arnold(InpD, ConS, Syn)


N1=length(InpD);
N2=length(ConS);


%%%%%%%%%%%%%%%
Inpval=unique(InpD);
Conval=unique(ConS);

bin_cop=[0:0.008:0.081]% linspace(  min(Conval), max(Conval),15);
bin_cop_width=bin_cop(2)-bin_cop(1);
bin_freq= -3.5:0.3:3.5;%linspace(  min(Inpval), max(Inpval),30);
bin_freq_width=bin_freq(2)-bin_freq(1);

for ind=1:length(bin_cop)
    for ind2=1:length(bin_freq)
   t=  find((ConS > bin_cop(ind)-bin_cop_width & ConS < bin_cop(ind)+bin_cop_width ) &   (InpD > bin_freq(ind2)-bin_freq_width & InpD < bin_freq(ind2)+bin_freq_width)    );
  if t>0
   arno_map(ind,ind2)= nanmean(Syn(t));
  else
         arno_map(ind,ind2)= NaN;
  end
    end
end




