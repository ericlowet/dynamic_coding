function [ TE12, TE21, TE11] =TE_sp_comp(ph1,ph2, t_s) 

% 

%% padding
ph1= [ zeros(1,1000) ph1 zeros(1,1000)];
ph2= [ zeros(1,1000) ph2 zeros(1,1000)];


  xt=ph1(1000:end-1000)';
  yt=ph2(1000:end-1000)';
  ytt=ph2(1000+t_s:end-1000+t_s)';
  
  %%%%%%%%%%%%%
  joint_prob= zeros(2,2,2);
 for ind= 0:1
 sel=find(xt== ind);
 %edges{1} = [ 0 1];
 % edges{2} = [ 0 1];
 [N,C]= hist3([ytt(sel) , yt(sel)] ,[2 2 ]);
 joint_prob(:,:,ind+1)= N;
 end
 countprob=joint_prob;
 NN=length(xt);
  joint_prob= joint_prob./NN;
  countmarg = sum(countprob,3);
   countmarg2 = sum(countmarg,1)./NN;
  sumTE=0;  sumTE2=0;nn=0;clear allk
  sumTE3=0;   sumTE4=0;
  for x=1:2
  
      for y=1:2
            yy=  sum(countmarg(:,y));
                 th=  countmarg(x,y)/NN;
                 tx= countmarg2(y) ;
                     sumTE3= sumTE3       - th* log(th);  
                        sumTE4= sumTE4         +( th* log(tx));  
          for z=1:2
              nn=nn+1;
               
              if (joint_prob(x,y,z) >0)
                       tn=  countmarg(x,y) /yy;
                      tk= countprob(x,y,z)/sum(countprob(:,y,z));
      
        
   sumTE= sumTE         + joint_prob(x,y,z)* log( tk);  
   sumTE2= sumTE2       -  joint_prob(x,y,z)* log(tn); 
     
              end
          end
      end
  end
  TE12=sumTE +   sumTE2;
   TE11=   sumTE4+sumTE3 ;
  
  %%%%%%%%%%%%%
  %%%%%%%%%%%%%
  %%%%%%%%%%%
  
  
% xt=ph2(1:end-t_s)';
%   yt=ph1(1:end-t_s)';
%   ytt=ph1(1+t_s:end)';
   xt=ph2(1000:end-1000)';
  yt=ph1(1000:end-1000)';
  ytt=ph1(1000+t_s:end-1000+t_s)';
  
  joint_prob= zeros(2,2,2);
 for ind= 0:1
 sel=find(xt== ind);
 %edges{1} = [ 0 1];
 % edges{2} = [ 0 1];
 [N,C]= hist3([ytt(sel) , yt(sel)] ,[2 2 ]);
 joint_prob(:,:,ind+1)= N;
 end
 countprob=joint_prob;
  joint_prob= joint_prob./length(xt);
  countmarg = sum(countprob,3);
  sumTE=0;  sumTE2=0;nn=0;clear allk
  for x=1:2
      for y=1:2
              yy=  sum(countmarg(:,y));
          for z=1:2
              nn=nn+1;

              if (joint_prob(x,y,z) >0)
                      tk= countprob(x,y,z)/sum(countprob(:,y,z));
                tn=  countmarg(x,y)/yy;
                  
   sumTE= sumTE          + joint_prob(x,y,z)* log( tk);  
   
  % allk(nn) = tk;
   sumTE2= sumTE2       -   joint_prob(x,y,z)* log( tn); 
              end
          end
      end
  end
  TE21=sumTE +   sumTE2;
   TE22=sumTE2;
  
  %%%%%%%%%%%%%
  %%%%%%%%%%%%%
  %%%%%%%%%%%
  
  