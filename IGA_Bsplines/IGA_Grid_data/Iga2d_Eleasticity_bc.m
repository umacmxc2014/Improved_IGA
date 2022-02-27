function [A,rhs]=Iga2d_Eleasticity_bc(A,rhs,m,n,uNoEs,vNoEs,g,Coordinate,knotSpanIndex,Qw,ConPts,U,V,wights,Ubar,pu,Vbar,pv,up,down,left,right)

udof=2*m;
temp=zeros(1,udof);
temp(1:2:udof-1)=1;

Node_up=m*(n-1)+(1:m);%---- v=1对应的边；
tmp=zeros(1,udof);
tmp(1:2:udof-1)=Node_up;
tmp(2:2:udof)=Node_up;
Node_up=2*tmp-temp;


Node_down=1:m;%------ v=0对应的边；
% ------ tmp=zeros(1,udof);
tmp(1:2:udof-1)=Node_down;
tmp(2:2:udof)=Node_down;
Node_down=2*tmp-temp;



vdof=2*n;
temp=zeros(1,vdof);
temp(1:2:end-1)=1;
Node_left=1+(0:n-1)*m;
tmp=zeros(1,vdof);
tmp(1:2:vdof-1)=Node_left;
tmp(2:2:vdof)=Node_left;
Node_left=2*tmp-temp;%-------- u=0对应的节点编号；


Node_right=(1:n)*m;
tmp(1:2:vdof-1)=Node_right;
tmp(2:2:vdof)=Node_right;
Node_right=2*tmp-temp;%-------  u=1对应的节点编号；

if nargin==4
be=[Node_up,Node_down,Node_left,Node_right];
be=unique(be);
A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
A(be,be)=eye(length(be));
return
end

%--------------------下面来求与区域边界相交的单元编号--------
% -----  e=i1+(j1-1)*uNoEs;%====参数区域里的第e号单元；

Element_up=(1:uNoEs)+(vNoEs-1)*uNoEs;%-----与区域的上边界相邻的单元编号；
Element_down=1:uNoEs;%-------与区域下边界相邻的单元编号；
Element_left=1+(0:vNoEs-1)*uNoEs;%----与区域左边界相交的单元编号；
Element_right=(1:vNoEs)*uNoEs;%-------与区域右边界相邻的单元编号；

[mu,nv,ndim]=size(ConPts);

udof=2*(pu+1);
temp=zeros(1,udof);
temp(1:2:udof-1)=1;

np=4;

if up==1 %------ v=1对应的边；
	for i=1:uNoEs
	ConPts_u=ConPts(:,nv,:);%--------v=1;
	ConPts_u=reshape(ConPts_u,mu,ndim);ConPts_u=ConPts_u';
	wu=wights(:,nv)';  Qwu=Qw(:,end)';
	e=Element_up(i);
	ue=Coordinate(e,[1,2]);
	uspan=knotSpanIndex(e,1);
	row=(uspan-pu:uspan)+(n-1)*m;
	tmp=zeros(1,udof);
	tmp(1:2:udof-1)=row;
	tmp(2:2:udof)=row;
	row=2*tmp-temp;
	
	x1={ConPts_u,U,wu,pu,Qwu,Ubar,uspan,g};x2={np,ue(1),ue(2)};
	Fe=Gauss_1d(@(u)temp_u_rhs(u,x1{:}),x2{:});%----- s=Gauss_1d(f,np,a,b)
		rhs(row)=rhs(row)+Fe;
		end
		be=[Node_down,Node_left,Node_right];be=unique(be);
		A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
		A(be,be)=eye(length(be));
		end
		
	
		
		
if down==1%------- v=0;
	for i=1:uNoEs
	
	ConPts_u=ConPts(:,1,:);%--------v=0;
	ConPts_u=reshape(ConPts_u,mu,ndim);ConPts_u=ConPts_u';
	wu=wights(:,1)';Qwu=Qw(:,1)';
		e=Element_down(i);
	ue=Coordinate(e,[1,2]);
	uspan=knotSpanInspan(e,1);
	
	row=uspan-pu:uspan;
	tmp=zeros(1,udof);
	tmp(1:2:udof-1)=row;tmp(2:2:udof)=row;
	row=2*tmp-temp;
	x1={ConPts_u,U,wu,pu,Qwu,Ubar,uspan,g};x2={np,ue(1),ue(2)};
	Fe=Gauss_1d(@(u)temp_u_rhs(u,x1{:}),x2{:});%----- s=Gauss_1d(f,np,a,b)
		rhs(row)=rhs(row)+Fe;
		end
		be=[Node_up,Node_left,Node_right];be=unique(be);
		A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
		A(be,be)=eye(length(be));
		end
		

	
vdof=2*(pv+1);
temp=zeros(1,vdof);
temp(1:2:vdof-1)=1;
		
		
if left==1%------- u=0;
	for i=1:vNoEs
		ConPts_v=ConPts(1,:,:);%--------u=0;
	ConPts_v=reshape(ConPts_v,nv,ndim);ConPts_v=ConPts_v';
	wv=wights(1,:);Qwv=Qw(1,:);
		e=Element_left(i);
	ve=Coordinate(e,[3,4]);
	vspan=knotSpanInspan(e,2);
	
	row=1+((vspan-pv:vspan)-1)*m;
	tmp=zeros(1,vdof);
	tmp(1:2:vdof-1)=row;tmp(2:2:vdof)=row;
	row=2*tmp-temp;
	x1={ConPts_v,V,wv,pv,Qwv,Vbar,vspan,g};x2={np,ve(1),ve(2)};
	Fe=Gauss_1d(@(v)temp_v_rhs(v,x1{:}),x2{:});%----- s=Gauss_1d(f,np,a,b)
		rhs(row)=rhs(row)+Fe;
		end
		be=[Node_up,Node_down,Node_right];be=unique(be);
		A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
		A(be,be)=eye(length(be));
		end

		
	if right==1%------- u=1;
	for i=1:vNoEs
		ConPts_v=ConPts(mu,:,:);%--------u=0;
	ConPts_v=reshape(ConPts_v,nv,ndim);ConPts_v=ConPts_v';
	wv=wights(end,:);Qwv=Qw(end,:);
	e=Element_right(i);
	ve=Coordinate(e,[3,4]);
	vspan=knotSpanInspan(e,2);
	row=(vspan-pv:vspan)*m;
	tmp=zeros(1,vdof);
	tmp(1:2:vdof-1)=row;tmp(2:2:vdof)=row;
	row=2*tmp-temp;
	x1={ConPts_v,V,wv,pv,Qwv,Vbar,vspan,g};x2={np,ve(1),ve(2)};
	Fe=Gauss_1d(@(v)temp_v_rhs(v,x1{:}),x2{:});%----- s=Gauss_1d(f,np,a,b)
		rhs(row)=rhs(row)+Fe;
		end
		be=[Node_up,Node_down,Node_left];be=unique(be);
		A(be,:)=zeros;A(:,be)=zeros;rhs(be)=zeros;
		A(be,be)=eye(length(be));
    end
end


	
	
	
	function s=temp_u_rhs(u,ConPts,U,wights,p,Qw,Ubar,uspan,g)
	[W,~,C,DC]=NurbsCurve(ConPts,U,wights,p,u);
	 Nu=bsplinebasis(Ubar,p,u);
	udof=2*(p+1);
	R=zeros(udof,2);
	R_u=Qw(uspan-p:uspan)'.*Nu/W;
	R(1:2:udof-1,1)=R_u;
	R(2:2:udof,2)=R_u;
	s=R*g(C(1))*norm(DC);
	end
	
	
	
	function s=temp_v_rhs(v,ConPts,V,wights,q,Qw,Vbar,vspan,g)
	[W,~,C,DC]=NurbsCurve(ConPts,V,wights,q,v);
	Nv=bsplinebasis(Vbar,q,v);
	vdof=2*(q+1);
	R=zeros(vdof,2);
	R_v=Qw(vspan-q:vspan)'.*Nv/W;
	R(1:2:vdof-1,1)=R_v;
	R(2:2:vdof,2)=R_v;
	s=R*g(C(2))*norm(DC);
	end
	