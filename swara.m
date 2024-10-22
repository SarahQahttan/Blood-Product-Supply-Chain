function [CC results]=swara(C)
I=[0.95 0.05 0.05
    0.75 0.25 0.25
    0.50 0.50 0.50 
    0.25 0.75 0.75
    0.05 0.95 0.95];
[n m]=size(C);
[CI]=lingu_Mat(C,I);
for i=1:m
ag1(i)=agg({CI{:,i}});
S(i)=score(ag1{i});
end
CC=[cell2mat(CI);cell2mat(ag1)];
r=rankWithDuplicates(S);
[a b]=sort(S,'descend');
for i=1:m
  if i==1
      sig(i)=1;
      k(i)=1;
      rho(i)=1;
  else
        sig(i)=a(i)/a(i-1);
    k(i)=sig(i)+1;
    rho(i)=rho(i-1)/k(i);
  end
end
w=rho./sum(rho);
results=[S;sig(r);k(r);rho(r);w(r)];