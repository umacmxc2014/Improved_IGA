function s=dgquad(f,a1,b1,a2,b2)
gp=[0.9491079123 0.7415311856 0.4058451514]; gp=[-gp gp 0];   %——————积分节点；
gw=[0.1294849662 0.2797053915 0.3818300505]; gw=[gw gw 0.4179591837];%————————高斯积分对应的权系数；
h1=(b1-a1)/2;m1=(a1+b1)/2;
h2=(b2-a2)/2;m2=(a2+b2)/2;
Fx=@(X)h1*X+m1;
Fy=@(Y)h2*Y+m2;
x=Fx(gp);
y=Fy(gp);
ngp=length(gw);
s=f(0,0);
[m,n]=size(s);
s=zeros(m,n);
for i=1:ngp
	for j=1:ngp
		s=s+gw(i)*gw(j)*f(x(i),y(j));
	end
	end
	s=s*h1*h2;