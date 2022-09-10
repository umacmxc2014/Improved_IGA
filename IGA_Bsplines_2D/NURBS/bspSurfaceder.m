function [DSu,DSv]=bspSurfaceder(P,U,p,u,V,q,v)
%============== 计算B样条曲面(x(u,v),y(u,v),z(u,v))=S(u,v)在点(u,v)处的关于u，v的一阶偏导数值；
%=========矩阵P为控制点组成的m*n*ndim阶矩阵；
%======== [m,n,ndim]=size(P);
%-----------  p为u方向B样条基函数的次数；U为开型节点向量；共有m+p+1个点构成；基函数个数为m;
%%-----------  q为v方向B样条基函数的次数；V为开型节点向量；共有n+q+1个点构成；基函数个数为n;

uspan=findspan(U,p,u);%----- u方向上节点张成区间的指标；
vspan=findspan(V,q,v);
ndim=size(P,3);%======= 控制点的所在的空间维数，为2或者3；
temp=P(uspan-p:uspan,vspan-q:vspan,:);%====== 由基函数不为0时，对应的控制点构造计算关于u,v的一阶偏导的控制点矩阵；
ConPtsUbar=zeros(p,q+1,ndim);%==存储计算u方向的一阶偏导时需要的控制点矩阵，这时u方向只要p个控制点；v方向需要q+1个；
ConPtsVbar=zeros(p+1,q,ndim);%==存储计算v方向的一阶偏导时需要的控制点矩阵，这时u方向只要p+1个控制点；v方向需要q个；
tempU=(U(uspan+1:uspan+p)-U(uspan-p+1:uspan))';
tempU=tempU*ones(1,q+1);%======tempU为p*(q+1)阶矩阵，存储u方向计算关于u的一阶偏导时需要的新控制点时的对应节点差分；
tempV=V(vspan+1:vspan+q)-V(vspan-q+1:vspan);
tempV=ones(p+1,1)*tempV;%======= tempV为(p+1)*q阶矩阵，存储v方向上的计算关于v的一阶偏导时需要的新控制点节点差分；
for  i=1:ndim	
    ConPtsUbar(:,:,i)=p*(temp(2:end,:,i)-temp(1:end-1,:,i))./tempU;
    ConPtsVbar(:,:,i)=q*(temp(:,2:end,i)-temp(:,1:end-1,i))./tempV;
end

%% 计算S(u,v)关于u方向的一阶偏导数；
Ubar=U;
% Ubar(1)=[];Ubar(end)=[];
 Ubar([1,end])=[];
Nu=bsplinebasis(Ubar,p-1,u);%===== 计算u方向上新的不为0的p个B样条基函数，由于对u求了一阶导数，因此是p-1次；
Nv=bsplinebasis(V,q,v);
DSu=zeros(ndim,1);%====DSu存储S(u,v)关于u的一阶偏导数值；
for i=1:ndim
    temp=reshape(ConPtsUbar(:,:,i),p,q+1);
   DSu(i)=Nu'*temp*Nv;
	end

% for i=1:p
%	for j=1:q+1
%	temp=reshape(ConPtsUbar(i,j,:),ndim,1);
% DSu=DSu+temp*Nu(i)*Nv(j);
% end
% end

%% 计算S(u,v)关于v方向的一阶偏导数；
Vbar=V;
% Vbar(1)=[];Vbar(end)=[];
Vbar([1,end])=[];
Nu=bsplinebasis(U,p,u);
Nv=bsplinebasis(Vbar,q-1,v);%===== 计算v方向上新的不为0的q个B样条基函数，由于对v求了一阶导数，因此是q-1次；
DSv=zeros(ndim,1);%===== DSv存储S(u,v)关于v的一阶偏导数值；
for i=1:ndim
    temp=reshape(ConPtsVbar(:,:,i),p+1,q);
   DSv(i)=Nu'*temp*Nv;
end

% for i=1:p+1
%	for j=1:q
%	temp=reshape(ConPtsVbar(i,j,:),ndim,1);
% DSv=DSv+temp*Nu(i)*Nv(j);
% end
% end
%% end
end

%%
%=============================Test======================
% U=[0 0 0 1/2 1 1 1];p=2;u=3/10;m=4;
% V=[0 0 0 1 1 1];q=2;v=6/10;n=3;
% P=zeros(m,n,3);
% P(:,:,1)=[0 3 6 9;0 3 6 9;0 3 6 9]';
% P(:,:,2)=[0 0 0 0;2 2 2 2;4 4 4 4]';
% P(:,:,3)=[0 3 3 0;2 5 5 2;0 3 3 0]';
% [DSu,DSv]=bspSurfaceder(P,U,p,u,V,q,v)

%  输出结果为：
%  输出结果为：
% DS =
%   8.400000000000000  -0.000000000000000
%                   0   4.000000000000000
%   4.800000000000001  -0.799999999999999

%=======================================================
