function [Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(P,w,U,p,V,q,t)
if t==0
    Q=P;wbar=w;Ubar=U;Vbar=V;
end

if t>=1
[m,n,ndim]=size(P);
UBreaks=unique(U);NoUBreks=length(UBreaks);
VBreaks=unique(V);NoVBreks=length(VBreaks);

     mu=m+(NoUBreks-1)*t;
	
     Q=zeros(mu,n,ndim);
    wbar=zeros(mu,n);
 for j=1:n
	 temp=reshape(P(:,j,:),m,ndim);temp=temp';
      [Ubar,temp,wU]=DegreeElevCurve(temp,w(:,j)',U,p,t);temp=temp';
      Q(:,j,:)=temp;
      wbar(:,j)=wU';
 end

P=Q;w=wbar;
nv=n+(NoVBreks-1)*t;
Q=zeros(mu,nv,ndim);
wbar=zeros(mu,nv);
 
 for i=1:mu
	 temp=reshape(P(i,:,:),n,ndim);temp=temp';
      [Vbar,temp,wV]=DegreeElevCurve(temp,w(i,:),V,q,t);temp=temp';
      Q(i,:,:)=temp;
      wbar(i,:)=wV;
 end
 
end
end

 