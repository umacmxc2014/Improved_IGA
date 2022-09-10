function Pw=WightedConPtsSurface(P,w)
%===求计算NUBRS曲面在齐次坐标表示时，需要用到的带权系数的控制点矩阵；控制点矩阵P为m*n*ndim型
%== m、n分别为为u、v方向上的基函数的个数；ndim为控制点所在的空间维数，为2或者3；
%==== 输出量为Pw，是一个带权系数的控制点矩阵，为m*n*(ndim+1)型的矩阵；

[m,n,ndim]=size(P);
Pw=zeros(m,n,ndim+1);
Pw(:,:,end)=w;
for i=1:ndim
	Pw(:,:,i)=P(:,:,i).*w;
end
end
