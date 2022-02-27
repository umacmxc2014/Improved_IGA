function [W,DW,C,DC]=NurbsCurve(P,U,w,p,u)
%======== ����Nurbs�����ڵ�u����һ�׵���ֵ��
%====== PΪ��ndim*n�׾��󹹳ɵĿ��Ƶ����;w ��Ȩϵ����������һ��1*n�׵���������pΪnurbs�������ߵĴ�����
%========= UΪ���ͽڵ���������һ��1*m �׵���������
% i=findspan(U,p,u);
[ndim,~]=size(P);
% n=length(w);����������Ϊn;
Pw=WightedConPtsCurve(P,w);
Cw=PointOnbspCurve(Pw,U,p,u);
W=Cw(end);%====== Nurbs�����еķ�ĸW;
C=Project(Cw);
Cwder=bspCurveder(Pw,U,p,u);
Ader=Cwder(1:ndim);
DW=Cwder(end);
DC=(Ader-C*DW)/W;
end

%==================Test====================================================
%===== The Nurbs Book P126 Ex4.2; 
% U=[0 0 0 1 1 1];
% p=2;w=[1 1 2];
% P=[1 1 0;
%    0 1 1];u=0;
% derC=DerNurbsCurve(P,U,w,p,u)
% ������Ϊ��
% derC =
%  0
%  2
%  ��P126ҳ(4.9)��p*w(2)/(U(p+2)*w(1))*(P(:,2)-P(:,1))һ�£�
%==========================================================================