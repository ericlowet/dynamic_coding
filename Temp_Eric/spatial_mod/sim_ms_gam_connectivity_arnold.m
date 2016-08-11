% [ extent, probability, strength];
EE= [ 4 , 0.25,  0.1   ];
II=[  2 , 0.2,  0.4    ];
IE=[  2 , 0.2,  0.4   ];
EI=[  3 , 0.3,  0.9    ];

EE2=[  6,   0.8,   0.0   ];

E2E2=[ 2 , 0.2, 0.1    ];
I2I2=[ 2 , 0.2, 0.4    ];
I2E2=[ 2 , 0.2, 0.6    ];
E2I2=[ 2 , 0.2, 0.1    ];


%%Connectivity

sig=[EE(1) II(1) IE(1) EI(1) EE2(1) E2E2(1) E2I2(1) I2I2(1) I2E2(1) ]*dR;

p=[EE(2) II(2) IE(2) EI(2) EE2(2) E2E2(2) E2I2(2) I2I2(2) I2E2(2)  ];

strength=[EE(3) II(3) IE(3) EI(3) EE2(3) E2E2(3) E2I2(3) I2I2(3) I2E2(3) ];

[S, connMat_sep, R_e, R_i ] = genConnections_ericJune29(numNeurLin, dR, sig, p, strength,Ne2,Ni2);
S3=S.*100;
S=single(S);
%figure,imagesc(S)
%  figure,imagesc(reshape(S(780,1:1600),40,40))
for ing=1:Ne2
    RF_con(:,:,ing)=reshape(S(Ne1+ing,1:Ne1),numNeurLin,numNeurLin);
end
RF_con=nanmean(RF_con,3);RF_con_norm=RF_con./max(RF_con(:));
%figure,imagesc(reshape(S(1500,1:1600),40,40))
% 
% figure,
% subplot(2,2,1)
% imagesc(connMat_sep{1});colorbar
% subplot(2,2,2)
% imagesc(connMat_sep{3});colorbar
% subplot(2,2,3)
% imagesc(connMat_sep{4});colorbar
% subplot(2,2,4)
% imagesc(connMat_sep{2});colorbar



