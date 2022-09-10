function S=PointOnNurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v)
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

