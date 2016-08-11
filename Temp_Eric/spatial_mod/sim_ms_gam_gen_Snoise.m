    
 %sigI=3*dR;
k=var_E1.*randn(numNeurLin, numNeurLin);
k=reshape(k,numNeurLin,numNeurLin);
% Rdum=min(R_e):dR:max(R_e);
% Rdum=Rdum-mean(Rdum);filtKern=exp(Rdum.^2/(2*sigI^2));
%  filtKern=bsxfun(@times,filtKern,filtKern');
%  filtKern=filtKern/sum(filtKern(:));k=(ifft2(fft2(k).*fft2(filtKern)));
NoiseE= [1.*k(:)',var_E2.*randn(1,Ne-Ne1)];