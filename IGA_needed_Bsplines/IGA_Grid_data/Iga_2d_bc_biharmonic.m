function [A,rhs]=Iga_2d_bc_biharmonic(A,rhs,m,n)

Node_down=1:m;%------ v=0¶ÔÓ¦µÄ±ß£»
Node_up=m*(n-1)+(1:m);%---- v=1¶ÔÓ¦µÄ±ß£»
Node_left=1+(0:n-1)*m;
Node_right=(1:n)*m;
be=[Node_up,Node_down,Node_left,Node_right];



be=unique(be);
A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
A(be,be)=eye(length(be));

down_2nd = m + (1:m); % 最靠近底边的的那一层内部节点。
up_2nd = m*(n-2)+ (1:m);
left_2nd = 2 + m*((1:n)-1);
right_2nd = m-1 + m*((1:n)-1);

be_2nd =[down_2nd,up_2nd,left_2nd,right_2nd];
be_2nd=unique(be_2nd);
A(be_2nd,:)=zeros; A(:,be_2nd)=zeros;rhs(be_2nd)=zeros;
A(be_2nd,be_2nd)=eye(length(be_2nd));

