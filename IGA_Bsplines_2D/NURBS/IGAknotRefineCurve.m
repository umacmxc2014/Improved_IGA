function  [Ubar,wbar]=IGAknotRefineCurve(U,w,p,Refinement)

%===============������Ϣ��==================================================
%----����p��NUrbs���ߵĽڵ�����U����һЩ�½ڵ�X���ڵ�����UΪ1*m�͵���������
%---X��Ԫ�ظ���Ϊr�����մ�С�������У�����X��Ԫ�ض��ǽڵ�����U�е��ڵ㣻
%--- PΪndim*n�Ϳ��Ƶ����nΪNurbs�������ڵĿռ�ά����nΪ�����������������Ƶ������
%---- wΪ��ʼNurbs���ߵĻ�������Ӧ��Ȩϵ������һ��1*n�͵���������
%----- �����½ڵ�X��Ҫ��Nurbs���ߵ���״�Ͳ������ֲ��䣻
%==========================================================================

%========================== �����Ϣ=======================================
%----- UbarΪ����X��õ����µĽڵ���������һ��1*(m+r)����������
%---- wbarΪ�µ�Ȩϵ����������һ��1*(n+r)�͵���������
%----- QΪ�µĿ��Ƶ������һ��ndim*(n+r)�͵ľ���
%------ ʵ���ϣ���IGA�У�ֻ��Ҫ����Ϣ��Ubar �� wbar;
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