function [rank,y,n,np,nn,S,R,Q]=grey_vikor(C,w,MM,rho,nu)
[n1,mm]=size(C);
y=C./repmat(sqrt(sum(C.^2)),n1,1);
y0=max(y);
for i=1:mm
    if MM(i)==1
        y0(i)=max(y(:,i));
    else
        y0(i)=min(y(:,i));
    end
end
n=(min(min(abs(repmat(y0,n1,1)-y)))+rho*max(max(abs(repmat(y0,n1,1)-y))))./(abs(repmat(y0,n1,1)-y)+rho*max(max(abs(repmat(y0,n1,1)-y+0.00000001))));
[np,nn]=deal(zeros(1,mm));
for i=1:mm
    np(i)=max(n(:,i));
    nn(i)=min(n(:,i));
end
S=sum(repmat(w,n1,1).*((repmat(np,n1,1)-n)./(repmat(np,n1,1)-repmat(nn,n1,1))),2);
R=max((repmat(w,n1,1).*((repmat(np,n1,1)-n)./(repmat(np,n1,1)-repmat(nn,n1,1))))')';
a1=nu*(S-max(S))./(min(S)-max(S));
a2=(1-nu)*(R-max(R))./(min(R)-max(R));
Q=a1+a2;
rank=rankWithDuplicates(Q)';
end