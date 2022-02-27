function [Qw,Ubar,Vbar,dof]=IGAknotRefineSurface(wights,knotU,pu,knotV,pv,Refinement)
Ubar=knotU;Vbar=knotV;
[~,n]=size(wights);
for i=1:Refinement
	UBreks=unique(Ubar);VBreks=unique(Vbar);
     Xu=(UBreks(1:end-1)+UBreks(2:end))/2;
     Xv=(VBreks(1:end-1)+VBreks(2:end))/2;
     temp=[Ubar,Xu];Ubar=temp;Ubar=sort(Ubar); % Ubar=[Ubar,Xu];Ubar=sort(Ubar);
     temp=[Vbar,Xv];Vbar=temp;Vbar=sort(Vbar);% Vbar=[Vbar,Xv]; Vbar=sort(Vbar);
end

mu=length(Ubar)-pu-1;nv=length(Vbar)-pv-1;
dof=mu*nv;
Qw=zeros(mu,n);
for j=1:n
	[Ubar,wU]=IGAknotRefineCurve(knotU,wights(:,j)',pu,Refinement);
     Qw(:,j)=wU';
end
wights=Qw;
Qw=zeros(mu,nv);

for i=1:mu
	[Vbar,wV]=IGAknotRefineCurve(knotV,wights(i,:),pv,Refinement);
     Qw(i,:)=wV;
end
