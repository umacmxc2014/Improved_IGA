function [S,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v)
%=========计算Nurbs曲面S(u,v)在参数点(u,v)处关于u，v的一阶偏导数；
% S(u,v)=（x(u,v),y(u,v),z(u,v)）(三维空间的曲面)；
%--- p为u方向的B样条基函数的次数；q为v方向的B样条基函数的次数；
% m=length(U)-p-1;%==== x方向的节点向量U的基函数个数为m;
% n=length(V)-q-1;%====== y方向的节点向量的基函数个数为n ；
%========= 矩阵P为控制点构成的矩阵，是 m * n *ndim型的控制网； w为m*n阶权系数矩阵；

%---- 输出的DS为ndim*2型的矩阵；第一列为NURBS曲面S(u,v)关于u方向的一阶偏导数，
%---    --- 而第二列为关于v的一阶偏导数值，即每个坐标分量关于u，v的偏导数；

[~,~,ndim]=size(ConPts);
Pw=WightedConPtsSurface(ConPts,wights);
Sw=PointOnBspSurface(Pw,knotU,pu,u,knotV,pv,v);
S=Project(Sw);
W=Sw(end);%==== Nurbs曲面中的分母W(u,v);
[DSwu,DSwv]=bspSurfaceder(Pw,knotU,pu,u,knotV,pv,v);%=====计算带权系数的控制点构成的多项式曲面关于u,v的一阶偏导数DSwu 、DSwv ;
DAu=DSwu(1:ndim);DAv=DSwv(1:ndim);DWu=DSwu(end);DWv=DSwv(end);
DSu=(DAu-DWu*S)/W;% Nurbs曲面S(u,v)关于u的一阶偏导数值；
DSv=(DAv-DWv*S)/W;% Nurbs曲面S(u,v)关于v的一阶偏导数值；
DF=[DSu,DSv];
 