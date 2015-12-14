function [prob,n1,n2]=double_exp_rev(timed, C1,C2,sign)




 C=C1;
 clear n1
  nn=0;
for t=timed
        nn=nn+1;
        V=t;
        if  ((-V)-C(3))>0
n1(nn)=  C(1)+C(2)*(V-C(3))+(C(5).*(   exp( ((-V)-C(3))./C(4)  )   ));
        else
    n1(nn)= 0;   
        end
end
 C=C2;
 clear n2
 nn=0;
for t=timed
       V=t;
    nn=nn+1;
    if ((-V)-C(3))>0
n2(nn)=  C(1)+C(2)*(V-C(3))+(C(5).*(   exp( ((-V)-C(3))./C(4)  )   ));
    else
     n2(nn)= 0;
    end
end

prob=(n1-n2)';
prob=sign.*( prob./max(prob));