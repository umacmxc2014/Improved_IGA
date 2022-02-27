function S=PointOnBspSurface(P,U,p,u,V,q,v)
%%==============计算B样条曲面参数曲面（x(u,v),y(u,v),z(u,v)）= S(u,v)在参数点（u，v）处的函数值；
%--- p为u方向的B样条基函数的次数；q为v方向的B样条基函数的次数；
% m=length(U)-p-1;%==== x方向的节点向量U的基函数个数为m;
% n=length(V)-q-1;%====== y方向的节点向量的基函数个数为n ；
%========= 矩阵P为控制点构成的矩阵，是 m * n *ndim型的控制网； 
%%
ndim=size(P,3);
uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
S=zeros(ndim,1);

for i=1:ndim
  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
 S(i)=Nu'*temp*Nv;
end
%==其他处理方法可见PointOnBspSurface1与PointOnBspSurface2，方法是一致的，只是在处理分量上，这里显得自然一些；
% S=zeros(1,1,ndim);
% for i=0:p
% for j=0:q
%    row=uspan-p+i;column=vspan-q+j;
%    S=S + P(row,column,:)*Nu(i+1)*Nv(j+1);
% end
% end
% S=reshape(S,ndim,1);



%%
%%==========================Test==================
% 见于 The NuRBS Book Page 116
% U=[0 0 0 1/2 1 1 1];p=2;u=3/10;m=4;
% V=[0 0 0 1 1 1];q=2;v=6/10;n=3;
% P=zeros(m,n,3);
% P(:,:,1)=[0 3 6 9;0 3 6 9;0 3 6 9]';
% P(:,:,2)=[0 0 0 0;2 2 2 2;4 4 4 4]';
% P(:,:,3)=[0 3 3 0;2 5 5 2;0 3 3 0]';
% S=PointOnBspSurface(P,U,p,u,V,q,v)
% ========= 输出结果为：
% S=
 %  3.060000000000000
 %  2.400000000000000
 %  3.480000000000000
 %%===========================================
 %% function S=PointOnBspSurface1(P,U,p,u,V,q,v)

% 若P=reshape(P,m*n,ndim)，则可以用下面的算法实现；
%% 下面是当把控制点存入（m*n）*ndim型的矩阵中，行数为二元B样条曲面的基函数的的个数，即为m*n个，ndim为空间的维数；
% m=length(U)-p-1;%------ 由节点向量U中构造出来的B样条基函数的个数；
% n=length(V)-q-1;%----- 由节点向量V构造出的B样条基函数的个数；则二元基函数的总个数为m*n个；
%% ndim=size(P,2);%  ------- 控制点所在空间的维数；
%% uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u); 
%% vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
%% S=zeros(1,ndim);
%% for j=vspan-q:vspan
%% index=(j-1)*m;  
%% for i=uspan-p:uspan
%% S=S+P(index+i,:)*Nu(i-uspan+p+1)*Nv(j-vspan+q+1);
%% end
%% end
%% S=S';
%% end
 
 %% Test=====================
 % The Nurbs Book P116 3.8.
 % P=[0,0,0;3,0,3;6,0,3;9,0,0;0,2,2;3,2,5;6,2,5;9,2,2;0,4,0;3,4,3;6,4,3;9,4,0];
% U=[0 0 0 1/2 1 1 1];V=[0 0 0  1 1 1];
% p=2;q=2;
% u=3/10;v=6/10;
% S=PointOnBspSurface1(P,U,p,u,V,q,v)
% S =
%   3.060000000000000
%   2.400000000000000
%   3.480000000000000
 %%==========================
 
 %% 根据对矩阵分量处理的不同，下面给出另外一个小程序；
%% function S=PointOnBspSurface2(P,U,p,u,V,q,v)
%%==============
%--- p为u方向的B样条基函数的次数；q为v方向的B样条基函数的次数；
% m=length(U)-p-1;%==== x方向的节点向量U的基函数个数为m;
% n=length(V)-q-1;%====== y方向的节点向量的基函数个数为n ；
%========= 矩阵P为控制点构成的矩阵，是 m * n *ndim型的控制网； 

%% ndim=size(P,3);
%% uspan=findspan(U,p,u);Nu=bsplinebasis(U,p,u);
%% vspan=findspan(V,q,v);Nv=bsplinebasis(V,q,v);
%% S=zeros(ndim,1);
%%  for i=1:ndim
%%  temp=reshape(P(uspan-p:uspan,vspan-q:vspan,i),p+1,q+1);
%% S(i)=Nu'*temp*Nv;
%%  end
%%  end

