function [S,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v)
%=========����Nurbs����S(u,v)�ڲ�����(u,v)������u��v��һ��ƫ������
% S(u,v)=��x(u,v),y(u,v),z(u,v)��(��ά�ռ������)��
%--- pΪu�����B�����������Ĵ�����qΪv�����B�����������Ĵ�����
% m=length(U)-p-1;%==== x����Ľڵ�����U�Ļ���������Ϊm;
% n=length(V)-q-1;%====== y����Ľڵ������Ļ���������Ϊn ��
%========= ����PΪ���Ƶ㹹�ɵľ����� m * n *ndim�͵Ŀ������� wΪm*n��Ȩϵ������

%---- �����DSΪndim*2�͵ľ��󣻵�һ��ΪNURBS����S(u,v)����u�����һ��ƫ������
%---    --- ���ڶ���Ϊ����v��һ��ƫ����ֵ����ÿ�������������u��v��ƫ������

[~,~,ndim]=size(ConPts);
Pw=WightedConPtsSurface(ConPts,wights);
Sw=PointOnBspSurface(Pw,knotU,pu,u,knotV,pv,v);
S=Project(Sw);
W=Sw(end);%==== Nurbs�����еķ�ĸW(u,v);
[DSwu,DSwv]=bspSurfaceder(Pw,knotU,pu,u,knotV,pv,v);%=====�����Ȩϵ���Ŀ��Ƶ㹹�ɵĶ���ʽ�������u,v��һ��ƫ����DSwu ��DSwv ;
DAu=DSwu(1:ndim);DAv=DSwv(1:ndim);DWu=DSwu(end);DWv=DSwv(end);
DSu=(DAu-DWu*S)/W;% Nurbs����S(u,v)����u��һ��ƫ����ֵ��
DSv=(DAv-DWv*S)/W;% Nurbs����S(u,v)����v��һ��ƫ����ֵ��
DF=[DSu,DSv];
 