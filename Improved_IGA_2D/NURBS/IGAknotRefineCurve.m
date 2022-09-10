function  [Ubar,wbar]=IGAknotRefineCurve(U,w,p,Refinement)

%===============输入信息：==================================================
%----对于p次NUrbs曲线的节点向量U插入一些新节点X；节点向量U为1*m型的行向量；
%---X的元素个数为r，按照从小到大排列，而且X的元素都是节点向量U中的内点；
%--- P为ndim*n型控制点矩阵，n为Nurbs曲线所在的空间维数，n为基函数个数，即控制点个数；
%---- w为初始Nurbs曲线的基函数对应的权系数，是一个1*n型的行向量；
%----- 插入新节点X后，要求Nurbs曲线的形状和参数保持不变；
%==========================================================================

%========================== 输出信息=======================================
%----- Ubar为插入X后得到的新的节点向量，是一个1*(m+r)型行向量；
%---- wbar为新的权系数向量，是一个1*(n+r)型的行向量；
%----- Q为新的控制点矩阵，是一个ndim*(n+r)型的矩阵；
%------ 实际上，在IGA中，只需要的信息是Ubar 和 wbar;
%==========================================================================

if Refinement==0
	Ubar=U;wbar=w;
else
	for i=1:Refinement
		UBreks=unique(U);
	    n=length(w);
	    X=(UBreks(2:end)+UBreks(1:end-1))./2;
        
	    Ubar=[U,X];
	    Ubar=sort(Ubar);
        r=length(X);
wbar=zeros(1,n+r);

 a=findspan(U,p,X(1)); 
 b=findspan(U,p,X(end));
 b=b+1;
 wbar(1:a-p)=w(1:a-p);
 wbar(b+r-1:n+r)=w(b-1:n);

 i=b+p-1;k=b+p+r-1;
 
 for j=r:-1:1
	 while (X(j)<=U(i))&&(i>a)
		 wbar(k-p-1)=w(i-p-1);
	 k=k-1;i=i-1;
	 end
	 wbar(k-p-1)=wbar(k-p);
     
	 for l=1:p
		 ind=k-p+l;
	 alpha=Ubar(k+l)-X(j);
	 if(abs(alpha)==0)
		 wbar(ind-1)=wbar(ind);
     else
         alpha=alpha/(Ubar(k+l)-U(i-p+l));
	  wbar(ind-1)=alpha*wbar(ind-1)+(1-alpha)*wbar(ind);
	 end
     end
     
	 k=k-1;
 end
U=Ubar;w=wbar;
	
end
 
end