function Cder=bspCurveder(P,U,p,u)
%============== 计算B样条曲线C(u)在点u处的关于u的导数；
%=========矩阵P为控制点组成的ndim*n阶矩阵；
%======== [ndim,n]=size(P);
%-----------  p为B样条基函数的次数；U为开型节点向量；共有n+p+1个点构成；基函数个数为n;
uspan=findspan(U,p,u);%----- 节点张成区间的指标；
ndim=size(P,1);%======= 控制点的空间维数，为2或者3；
temp=P(:,uspan-p:uspan);
Q=zeros(ndim,p);%==存储计算关于u的一阶导数时需要的控制点矩阵，这时一共有ndim*p，即p个控制点；

%% 计算u方向的一阶导数；
for  i=1:ndim
	Q(i,:)=p*(temp(i,2:end)-temp(i,1:end-1))./(U(uspan+1:uspan+p)-U(uspan-p+1:uspan));
end
Ubar=U;
Ubar([1,end])=[];%-----------即 Ubar(1)=[];Ubar(end)=[];
Nu=bsplinebasis(Ubar,p-1,u);%===== 计算u方向上新的不为0的p个B样条基函数，由于对u求了一阶导数，因此是p-1次；
Cder=Q*Nu;
end

%% end
%------------------Test----------------
% P=[0,1,1,2 3 4 5 6;
%    1,1,0 2 3 4 5 6];
% p=2;U=[0 0 0 1 2 3 4 4 5 5 5];u=5/2;
%   Cder=bspCurveder(P,U,p,u)
% 输出结果为：
% Cder =
%   1.000000000000000
%   1.500000000000000
%----------------------------------------