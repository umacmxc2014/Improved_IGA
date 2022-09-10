function [Ubar,Q,wbar]=DegreeElevCurve(P,w,U,p,t)
[ndim,n]=size(P);
Pw=WightedConPtsCurve(P,w);
UBreks=unique(U);NoBreks=length(UBreks);
n=n+(NoBreks-1)*t;%-----对Nurbs曲线升阶t次后的基函数个数；
ndim=ndim+1;%----带权系数矩阵的空间维数为ndim;
Qw=zeros(ndim,n);
bezalfs=zeros(p+t+1,p+1);
bpts=zeros(ndim,p+1);%----p次Bezier曲线的p+1个控制点；
ebpts=zeros(ndim,p+t+1);%---- p+t次Bezier曲线的p+t+1个控制点；
Nextbpts=zeros(ndim,p-1);
alfs=zeros(1,p-1);
temp=length(U);
m=temp+NoBreks*t;
Ubar=zeros(1,m);%----------升阶t次后的新的节点向量元素的个数为m;
m=temp-1;%-------由于在The Nurbs Book中节点向量的个数为m+1，因此这里要把m变为m-1;

ph=p+t;ph2=fix(ph/2);
bezalfs(1,1)=1.0;bezalfs(ph+1,p+1)=1.0;
for i=1:ph2
	inv=1/Bin(ph,i);
mpi=min(p,i);

for j=max(0,i-t):mpi
	bezalfs(i+1,j+1)=inv*Bin(p,j)*Bin(t,i-j);
end
end

for i=ph2+1:ph-1
	mpi=min(p,i);
for j=max(0,i-t):mpi
	bezalfs(i+1,j+1)=bezalfs(ph-i+1,p-j+1);
end
end

mh=ph;kind=ph+1;
r=-1;a=p;
b=p+1;cind=1;
ua=U(1);Qw(:,1)=Pw(:,1);
Ubar(1:ph+1)=ua*ones(1,ph+1);
% for i=0:ph
% Uh(i+1)=ua;
% end
bpts(:,1:p+1)=Pw(:,1:p+1);
% for i=0:p
%	bpts(:,i+1)=Pw(:,i+1);
% end

while(b<m)
	i=b;
while(b<m && U(b+1)==U(b+2)) 
    b=b+1;
end
mul=b-i+1;
mh=mh+mul+t;
ub=U(b+1);
oldr=r;r=p-mul;

if(oldr>0)
    lbz=fix((oldr+2)/2);
else
lbz=1;
end

if(r>0)
    rbz=ph-fix((r+1)/2);
else
	rbz=ph;
end

if(r>0) 
	numer=ub-ua;
for k=p:-1:(mul+1)
	alfs(k-mul)=numer/(U(a+k+1)-ua);
end

for j=1:r
	save=r-j;s=mul+j;
    
for k=p:-1:s
	bpts(:,k+1)=alfs(k-s+1)*bpts(:,k+1)+(1-alfs(k-s+1))*bpts(:,k);
end
Nextbpts(:,save+1)=bpts(:,p+1);
end
end

for i=lbz:ph
ebpts(:,i+1)=zeros;mpi=min(p,i);

for j=max(0,i-t):mpi
	ebpts(:,i+1)=ebpts(:,i+1)+bezalfs(i+1,j+1)*bpts(:,j+1);
end
end

if(oldr>1)
first=kind-2;last=kind;
den=ub-ua;
bet=fix((ub-Ubar(kind))/den);

for tr=1:(oldr-1)
	i=first;j=last;kj=j-kind+1;
    
while(j-i>tr)
    
	if(i<cind)
		alf=(ub-Ubar(i+1))/(ua-Ubar(i+1));
	Qw(:,i+1)=alf*Qw(:,i+1)+(1-alf)*Qw(:,i);
    end
    
	if(j>=lbz)
        
		if(j-tr<=kind-ph+oldr)
			gam=(ub-Ubar(j-tr+1))/den;
		ebpts(:,kj+1)=gam*ebpts(:,kj+1)+(1-gam)*ebpts(:,kj+2);
		else
		ebpts(:,kj+1)=bet*ebpts(:,kj+1)+(1-bet)*ebpts(:,kj+2);
        end
        
		end
		i=i+1;j=j-1;kj=kj-1;
		end
		first=first-1;last=last+1;
		end
end
        
		if a~=p
			for i=0:(ph-oldr-1)
				Ubar(kind+1)=ua;kind=kind+1;
            end
        end
            
			for j=lbz:rbz
				Qw(:,cind+1)=ebpts(:,j+1);cind=cind+1;
            end
            
			if(b<m)
				for j=0:(r-1)
					bpts(:,j+1)=Nextbpts(:,j+1);
				end
	 bpts(:,(r:p)+1)=Pw(:,b-p+r+1:b+1); 
            % for j=r:p
			%	bpts(:,j+1)=Pw(:,b-p+j+1);
            % end
				a=b;b=b+1;ua=ub;
				
				else
		Ubar(kind+1:kind+ph+1)=ub*ones(1,ph+1);
        	%		for i=0:ph
            % Uh(kind+i+1)=ub;
			%	end
            end
end
wbar=Qw(end,:);
Q=Qw(1:ndim-1,:);
for i=1:ndim-1
    Q(i,:)=Q(i,:)./wbar;
end
%================================Test=====================================
% P=[1 0;1 1;0 1]';p=2;w=[1 1 2];
% U=[0 0 0 1 1 1];t=2;u=1/4;
% C=PointOnNurbsCurve(P,U,w,p,u)
% C =
%   0.882352941176471
%   0.470588235294118
%---------------------升阶t后，求出新的控制点Q和权系数wbar,Ubar,带入计算-----
% [Ubar,Q,wbar]=DegreeElevCurve(P,w,U,p,t)
% C=PointOnNurbsCurve(Q,Uh,wbar,p+t,u)
% C =
%   0.882352941176471
%   0.470588235294118
%----------------计算结果一致----------------------------------------------
%==========================================================================
					