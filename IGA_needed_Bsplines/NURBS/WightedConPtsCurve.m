function Pw=WightedConPtsCurve(P,w)
%===求计算NUBRS曲线在齐次坐标表示时，需要用到的带权系数的控制点矩阵；控制点矩阵P为ndim*n型;
%== n为控制点的个数，即u方向上基函数的个数为n；
%=== 输出量为Pw是一个(ndim+1)*n型的带权系数的控制点矩阵；

[ndim,n]=size(P);
Pw=zeros(ndim+1,n);
Pw(end,:)=w;
for i=1:ndim
	Pw(i,:)=P(i,:).*w;
end
    
end