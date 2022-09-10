function  N=bsplinebasis(U,p,u)
% Compute the (p+1) p-th B-splines basis functions that are non-zero at the parameter value u.
N=zeros(p+1,1);%The column vector N stores the non-zero B-spline basis functions
N(1)=1.0;
i=findspan(U,p,u);% the knot span index of u, $u\in [U_{i}, U_{i+1})$.

for j=1:p
  saved=0;
for k=1:j
temp=N(k);
% N(k)  = saved + temp*(U(i+k)-u)/(U(i+k)-U(i+k-j));
% saved = temp*(u-U(i+k-j))/(U(i+k)-U(i+k-j));
left=U(i+k)-u;  right=u-U(i+k-j);
tmp=temp/(U(i+k)-U(i+k-j));
N(k)=saved+tmp*left;
saved=tmp*right;
end
N(j+1)=saved;
end
end

%======================================== Test=============================
%------------- The Nurbs----- P69,Ex2.3---------
% U=[0,0,0,1,2,3,4,4,5,5,5];
% p=2;u=5/2;
% N=bsplinebasis(U,p,u)'
% ==== 输出结果为========
% N =
%  0.125000000000000   0.750000000000000   0.125000000000000
%==========================================================================
