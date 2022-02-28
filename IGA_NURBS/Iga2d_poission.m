function [err,dof]=Iga2d_poission(ConPts,weights,knotU,pu,knotV,pv,Refinement,t, test_case)
% This numerical example is taken from Paper: 
%  Numerical Demonstration of Finite Element Convergence  forLagrange Elements in COMSOL Multiphysics
%  Section  4   Curved Boundary;
%  The physical domain is a unit disk ;

%=====================
% Input:
% ConPts are the control points;
% weights are the weights;
% knotU contains the knot vector in the u direction;
% pu is the degree of NURBS surface in the u direction;
% knotV contains the knot vector in the v direction;
% pv is the degree of NURBS surface in the v direction;
% Refinement mens the times of the h-refinement of both u and v direction;
% t denotes the order to be elevated with analogous to p-refin, i.e., the
% ultimate degree of NURBS  basis functions is (pu+t) in u direction, and
% (pv+t) in the v direction;

addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')



if(strcmp(test_case ,'rectangle'))

u_Exact=@(x,y)sin(pi*x)*(y.^2-1);% Exact solution of Poisson equation;
f=@(x,y) -(-pi*pi*sin(pi*x).*(y.^2-1)+2*sin(pi*x)) ;% The right hand side of Poisson equation;
u_d=@(x,y)sin(pi*x)*(y.^2-1);
u_grad=@(x,y)[ pi*(y^2-1)*cos(pi*x), 2*y*sin(pi*x) ];

end

if(strcmp(test_case,'quarter'))
       
u_Exact=@(x,y)  x.^2.*y.^2*sin(pi*(x.^2+y.^2-2));% Exact solution of Poisson equation;
f=@(x,y) -(  20*pi*x.^2.*y.^2.*cos(pi*(x.^2+y.^2-2)) + ...
      2*(x^2+y^2-2*pi*pi*x^4*y^2-2*pi*pi*x^2*y^4)*sin(pi*(x.^2+y.^2-2))  );% The right hand side of Poisson equation;
u_d=@(x,y) x.^2.*y.^2*sin(pi*(x.^2+y.^2-2));

u_grad=@(x,y) [ y^2*(2*x*sin(pi*(x.^2+y.^2-2)) +2*pi*x^3*cos(pi*(x.^2+y.^2-2))) , ...
               x^2*(2*y*sin(pi*(x.^2+y.^2-2)) +2*pi*y^3*cos(pi*(x.^2+y.^2-2))) ] ;
    
       
% u_Exact=@(x,y)  x.^2.*y.^2*(x.^2+y.^2-2);% Exact solution of Poisson equation;
% f=@(x,y) -(20*x^2*y^2 + 2*x^2*(x^2 + y^2 - 2) + 2*y^2*(x^2 + y^2 - 2));% The right hand side of Poisson equation;
% u_d=@(x,y) x.^2.*y.^2*(x.^2+y.^2-2);

% u_grad=@(x,y) [y^2*(2*x*(x^2+y^2-2) +2*x^3), ...
%               x^2*(2*y*(y^2+x^2-2) +2*y^3)] ;
    
end


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
 ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end

nurbsInfo=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);


Element=nurbsInfo.Element;
Coordinate=nurbsInfo.Coordinate;
knotSpanIndex=nurbsInfo.knotSpanIndex;
Ubar=nurbsInfo.Ubar;
Vbar=nurbsInfo.Vbar;
m=nurbsInfo.m;
n=nurbsInfo.n;
Qw=nurbsInfo.Qw;
NoEs=nurbsInfo.NoEs;
dof=nurbsInfo.dof;

Ubreaks=nurbsInfo.UBreaks;
Vbreaks=nurbsInfo.VBreaks;

uNoEs=nurbsInfo.uNoEs;
vNoEs=nurbsInfo.vNoEs;


[n_conpts_u,n_conpts_v,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

U_ConPts=zeros(DIM,n_conpts_u);

for i=1:DIM
  for j=1:n_conpts_u
    U_ConPts(i,j)=ConPts(j,1,i);
end
end

U_weights=weights(:,1)';

U_wbar=Qw(:,1);

 [err,u_d_h] = L2_project2dirichlet_bnd(U_ConPts,knotU,U_weights,pu, Ubar,U_wbar, Ubreaks, uNoEs, u_d);
 
 disp('L2 projection error is ')
 disp(err)

% [err,u_d_h] = L2_project2dirichlet_bnd(ConPts,weights,knotU,pu,knotV,pv, Element,Coordinate, knotSpanIndex, Qw, m,Ubar, Vbar, uNoEs, u_d);

% disp('The L2 projection error  is ')
% disp(u_d_h)

A=sparse(dof,dof);rhs=zeros(dof,1);
Eledof=(pu+1)*(pv+1);% The number of dof in "Element";

np=pu+1;% The number of Gauss quadrature points in  element;

if(np>=9)
    np=9;
end

disp('The degree of the  NURBS basis is ')
disp(pu)

disp('NoEs=')
disp(NoEs)

for e=1:NoEs % Loop for each element;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);
    uspan=knotSpanIndex(e,1);
   
    vspan=knotSpanIndex(e,2);
    m_row = Eledof;
    n_column = Eledof + 1; 
	x1={Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv,Eledof};
	x2={np,ue(1),ue(2),ve(1),ve(2),m_row,n_column};

    s=Gauss_2d(@(u,v)quad_Ae_Fe(u,v,x1{:},f),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)
  
	Ae=s(1:Eledof,1:Eledof);
    Fe=s(:,end);
	A(row,row)=A(row,row)+Ae;
    rhs(row)=rhs(row)+Fe;
end

% 首先处理下边界 (i.e., v=0) 是 非齐次 Dirichlet 边界。

ve=[Vbar(1), Vbar(pv+2)];  % The support of  R_{1,pv}(v) is [v_1,v_{1+pv+1}]; 

for e=1:uNoEs % Loop for the  elements along the boundary v=0;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    
    uspan=knotSpanIndex(e,1);
    m_row = Eledof;
    n_column = 1;
	x1={Qw,Ubar,uspan,Vbar,ConPts,weights,knotU,pu,knotV,pv,Eledof};
	x2={np,ue(1),ue(2),ve(1),ve(2),m_row,n_column};
    
    s=Gauss_2d(@(u,v)quad_Ae_Fe_u(u,v,x1{:},u_d_h ),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)

    rhs(row)=rhs(row)-s;
end



[A,rhs]=Iga_2d_bc(A,rhs,m,n);

% spy(A)

U=A\rhs;

% 这时候先把齐次解$u_0$算出来了。
% 现在把非齐次狄利克雷边界的解加上上面的非齐次解里。

U(1:m)=U(1:m)+u_d_h;

U=reshape(U,m,n);





%
err=[0,0];
	    
 for e=1:NoEs
    ue=Coordinate(e,1:2);ve=Coordinate(e,3:4);
    uspan=knotSpanIndex(e,1);vspan=knotSpanIndex(e,2);
    m_row = 1;
    n_column = 2;
	x1={U,Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv};
	x2={np,ue(1),ue(2),ve(1),ve(2),m_row,n_column};
    err=err+Gauss_2d(@(u,v)Err_L2(u,v,x1{:},u_Exact,u_grad),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)
 end
	err=sqrt(err);
end


function s=quad_Ae_Fe(u,v,Qw,Ubar,uspan,Vbar,vspan,ConPts,wights,knotU,pu,knotV,pv,Eledof,f)
% ----- Eledof=(pu+1)*(pv+1);
[F,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v);
Uders=bspbasisDers(Ubar,pu,u,1);Vders=bspbasisDers(Vbar,pv,v,1);
Nu=Uders(1,:)';DNu=Uders(2,:)';
Nv=Vders(1,:);DNv=Vders(2,:);
%--- i=uspan;j=vspan;
J=abs(det(DF));
temp=Qw(uspan-pu:uspan,vspan-pv:vspan);
R=temp.*(Nu*Nv)/W;R=reshape(R,Eledof,1);

W2=W^2;
s_Fe=f(F(1),F(2))*R*J;
DRu=temp.*((W*DNu-DWu*Nu)*Nv)/W2;DRu=reshape(DRu,Eledof,1);
DRv=temp.*(Nu*(DNv*W-DWv*Nv))/W2;DRv=reshape(DRv,Eledof,1);
DR=[DRu,DRv]/DF;

s_Ae=DR*DR'*J;
s=[s_Ae,s_Fe];
end




function s=quad_Ae_Fe_u(u,v,Qw,Ubar,uspan,Vbar,ConPts,wights,knotU,pu,knotV,pv,Eledof,u_d_h )
% ----- Eledof=(pu+1)*(pv+1);
[F,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v);
vspan=pv+1;
Uders=bspbasisDers(Ubar,pu,u,1);Vders=bspbasisDers(Vbar,pv,v,1);
Nu=Uders(1,:)';DNu=Uders(2,:)';
Nv=Vders(1,:);DNv=Vders(2,:);
%--- i=uspan;j=vspan;
J=abs(det(DF));
temp=Qw(uspan-pu:uspan,vspan-pv:vspan);


W2=W^2;

DRu=temp.*((W*DNu-DWu*Nu)*Nv)/W2;DRu=reshape(DRu,Eledof,1);
DRv=temp.*(Nu*(DNv*W-DWv*Nv))/W2;DRv=reshape(DRv,Eledof,1);
DR=[DRu,DRv]/DF;


DR_v0=DR(1:pu+1,:); % 这个是对应R_{j,1}(u,v)对应的基函数，对应于处理 v =0时下边界的非齐次狄利克雷边界条件时需要用到的NURBS基函数.

D_u_d_h_v0= u_d_h(uspan-pu:uspan)'*DR_v0*J; % 居然忘记乘以Jacobian 了。唉。

s= DR*D_u_d_h_v0';

end


function s=Err_L2(u,v,U,Qw,Ubar,uspan,Vbar,vspan,ConPts,wights,knotU,pu,knotV,pv,u_Exact, u_grad)


Uindex=uspan-pu:uspan;Vindex=vspan-pv:vspan;

[F,DF,W,DWu,DWv]=NurbsSurface(ConPts,wights,knotU,pu,u,knotV,pv,v);
Uders=bspbasisDers(Ubar,pu,u,1);Vders=bspbasisDers(Vbar,pv,v,1);
Nu=Uders(1,:)';DNu=Uders(2,:)';
Nv=Vders(1,:);DNv=Vders(2,:);
%--- i=uspan;j=vspan;
Eledof = (pu+1)*(pv+1);

J=abs(det(DF));
temp=Qw(uspan-pu:uspan,vspan-pv:vspan);
U_ij = U(Uindex,Vindex);
U_ij = U_ij(:);
U_ij = U_ij';

R=temp.*(Nu*Nv)/W;R=reshape(R,Eledof,1);
s_L2_err = (u_Exact(F(1),F(2)) - U_ij*R ).^2*J;

W2=W^2;

DRu=temp.*((W*DNu-DWu*Nu)*Nv)/W2;DRu=reshape(DRu,Eledof,1);
DRv=temp.*(Nu*(DNv*W-DWv*Nv))/W2;DRv=reshape(DRv,Eledof,1);
DR=[DRu,DRv]/DF;

s_H1_err = sum( (u_grad(F(1),F(2))  - U_ij*DR ).^2)*J;


s=[s_L2_err,s_H1_err];


end

% Example:
%=========================================================
%t=0;
% knotU=[0,0,0,1,1,1];
% knotV=knotU;
% a=sqrt(2)/2;
% wights=[1,a,1;a,1,a;1,a,1];
% ConPts=zeros(3,3,2);
% ConPts(:,:,1)=[-a,0,a; -2*a   0   2*a;-a, 0, a];
% ConPts(:,:,2)=[ a,2*a,  a; 0,0, 0;-a, -2*a, -a];
% pu=2;pv=2;
% Refinement=3;
% tic
% [err,dof]=Iga2d_poission(ConPts,wights,knotU,pu,knotV,pv,Refinement,t)
% toc
