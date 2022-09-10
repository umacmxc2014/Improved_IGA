function [A,rhs]=Iga_2d_bc(A,rhs,m,n)
Node_down=1:m;%------ v=0，即下边界。
Node_up=m*(n-1)+(1:m);%---- v=1,即上边界。
Node_left=1+(0:n-1)*m;
Node_right=(1:n)*m; 
 be=[Node_up,Node_down,Node_left,Node_right];
% be=[Node_left,Node_right];
be=unique(be);
A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
A(be,be)=eye(length(be));
