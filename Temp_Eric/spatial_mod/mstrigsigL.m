function  [ ms_sig] = mstrigsig(ms_dips,signals)

[n,m] =  size(signals);

ms_n=length(ms_dips);
win_l=450;
win_r=800;
n1=0;
for ind=1:ms_n
if ((ms_dips(ind)-win_l)>0)   &  ((ms_dips(ind)+win_r)< m)
n1=n1+1;
ms_sig(:,:,n1)=signals(:,ms_dips(ind)-win_l:ms_dips(ind)+win_r);
end
end










