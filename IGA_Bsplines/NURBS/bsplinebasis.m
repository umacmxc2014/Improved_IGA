function  N=bsplinebasis(U,p,u)
%=============计算p次B样条在节点向量U给定时，在参数点u处的p+1个不为0的基函数的函数值；即，N(i-p),...,N(i)
N=zeros(p+1,1);%--------用p+1维的列向量保存这p+1个不为0的基函数值；
N(1)=1.0;
i=findspan(U,p,u);%------求出点u所在的节点张成区间指标i，即u属于[U(i),U(i+1))

for j=1:p
	saved=0;
for k=1:j
temp=N(k);
% N(k)  = saved + temp*(U(i+k)-u)/(U(i+k)-U(i+k-j));
% saved = temp*(u-U(i+k-j))/(U(i+k)-U(i+k-j));
left=U(i+k)-u;  right=u-U(i+k-j);
tmp=temp/(U(i+k)-U(i+k-j));
N(k)=saved+tmp*left;
saved=tmp*right;
end
N(j+1)=saved;
end
end

%======================================== Test=============================
%------------- The Nurbs----- P69,Ex2.3---------
% U=[0,0,0,1,2,3,4,4,5,5,5];
% p=2;u=5/2;
% N=bsplinebasis(U,p,u)'
% ==== 输出结果为========
% N =
%  0.125000000000000   0.750000000000000   0.125000000000000
%==========================================================================