function [W,DW,C,DC]=NurbsCurve(P,U,w,p,u)
%======== 计算Nurbs曲线在点u处的一阶导数值；
%====== P为由ndim*n阶矩阵构成的控制点矩阵;w 是权系数向量，是一个1*n阶的行向量；p为nurbs样条曲线的次数；
%========= U为开型节点向量，是一个1*m 阶的行向量；
% i=findspan(U,p,u);
[ndim,~]=size(P);
% n=length(w);基函数个数为n;
Pw=WightedConPtsCurve(P,w);
Cw=PointOnbspCurve(Pw,U,p,u);
W=Cw(end);%====== Nurbs曲线中的分母W;
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
% 输出结果为：
% derC =
%  0
%  2
%  与P126页(4.9)：p*w(2)/(U(p+2)*w(1))*(P(:,2)-P(:,1))一致；
%==========================================================================