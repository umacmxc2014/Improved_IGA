function  [u_h_g1,u_h_g2,err] = biharmonic_L2_project2dirichlet_bnd(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,g1,g2)
% The dimension of ConPts is 

addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')

% g1=@(x,y) sin(pi*x+1)*sin(pi*y+1);% u = g1 on $\partial \Omega$.
% g2=@(x,y) sin(pi*x+1)*sin(pi*y+1); % grad u \cdot n =g2  on $\partial \Omega$.


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
 ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end


disp('pu=')
disp(pu)

nurbsInfo=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);



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
V_ConPts=zeros(DIM,n_conpts_v);

for i=1:DIM
  for j=1:n_conpts_u
    U_ConPts(i,j)=ConPts(j,1,i);
end
end

for i=1:DIM
  for j=1:n_conpts_v
    V_ConPts(i,j)=ConPts(j,2,i);
end
end

U_weights=weights(:,1)';

U_wbar=Qw(:,1);


 [err_g1, u_h_g1] = L2_project2dirichlet_bnd(U_ConPts,knotU,U_weights,pu, Ubar,U_wbar, Ubreaks, uNoEs, g1);
 
 ele_dof=pu+1;

coordinate=zeros(uNoEs,2);
element=zeros(uNoEs,ele_dof);
knotSpanIndex=zeros(uNoEs,1);

for i=1:uNoEs
coordinate(i,:)=[Ubreaks(i),Ubreaks(i+1)];
span=findspan(Ubar,pu,Ubreaks(i));
knotSpanIndex(i)=span;
element(i,:)=span-pu:span;
end

dofs=length(U_wbar);
A=sparse(dofs,dofs);
rhs=zeros(dofs,1);


np=pu+1;% 每一条边上的数值积分点的个数.
 
 for i=1:uNoEs
	  row=element(i,:);
     
    ue=coordinate(i,:);
    uspan=knotSpanIndex(i);
    m_row = ele_dof;
    n_column = ele_dof + 1;
    x1={np,ue(1),ue(2),m_row,n_column};
    s=Gauss_1d(@(u)quad_Ae_Fe_biharmonic(u,Ubar,Vbar, Qw,  ConPts,weights,knotU,knotV, pu,pv, uspan, g2, u_h_g1),x1{:});
     
    A(row,row)=A(row,row)+s(:,1:ele_dof);
    rhs(row)=rhs(row)+s(:,end);
end
 
 u_h_g2=A\rhs;
 
 
  err_g2=[0,0];

 for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 2;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_biharmonic(u,Ubar,Vbar, Qw,  ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, g2),x1{:});
   err_g2 = err_g2 + s;
 end
   
   err_g2=sqrt(err_g2);

   err=[err_g1,err_g2];

end
 
 
 


 
 function s=quad_Ae_Fe_biharmonic(u,Ubar,Vbar, Qw,  ConPts,weights,knotU,knotV, pu,pv, uspan, g2,u_h_g1)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。



[F,DF,W,DWu,DWv]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,0);




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



[W,DW,C,DC]=NurbsCurve(U_ConPts,knotU,U_weights,pu,u);



Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,0,1); %  这时候只需要 v =0 处的基函数及其导数值.

Nu=Uders(1,:)';
DNu=Uders(2,:)';

Nv=Vders(1,:);
DNv=Vders(2,:);



W2=W^2;
temp_2=Qw(uspan-pu:uspan,2);

DRu=temp_2.*((W*DNu-DWu*Nu)*Nv(2))/W2; 
DRv=temp_2.*(Nu*(DNv(2)*W-DWv*Nv(2)))/W2; 
DR=[DRu,DRv]/DF;


normal=zeros(2,1);
normal(1)= DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);

J=norm(DC,2);

DR_n=DR*normal;% 这时候是列向量.
DR_n=DR_n'; % 这时候的DR_n实际上是只对 构造g2的插值函数的逼近而已, 也就是 $\sum_{i=1}^{m} u_{i,2}R_{i,2}$的基函数,
% 即， $R_{i,2}$的梯度在 $v=0$时的法向分量.

temp=U_wbar(uspan-pu:uspan);
R=temp.*Nu/W;

Ae=R*DR_n*J;





% disp('J=')
% disp(J)




temp_1=Qw(uspan-pu:uspan,1); % 现在处理 u=u_g1时的L2投影的梯度点乘n与 v_h 的内积.
DRu=temp_1.*((W*DNu-DWu*Nu)*Nv(1))/W2; 
DRv=temp_1.*(Nu*(DNv(1)*W-DWv*Nv(1)))/W2; 
DR=[DRu,DRv]/DF;
DR_n=DR*normal;


D_u_g1_n= (u_h_g1(uspan-pu:uspan)' * DR_n);


Fe=(g2(C(1),C(2)) - D_u_g1_n)*R*J;

s=[Ae,Fe];

 end

 
 
 function s=quad_err_L2_biharmonic(u,Ubar,Vbar, Qw,  ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, g2)

[F,DF,W,DWu,DWv]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,0);

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



[W,DW,C,DC]=NurbsCurve(U_ConPts,knotU,U_weights,pu,u);



Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,0,1);

Nu=Uders(1,:)';
DNu=Uders(2,:)';

Nv=Vders(1,:);
DNv=Vders(2,:);



W2=W^2;
temp=Qw(uspan-pu:uspan,2);

DRu=temp.*((W*DNu-DWu*Nu)*Nv(2))/W2; 
DRv=temp.*(Nu*(DNv(2)*W-DWv*Nv(2)))/W2; 
DR=[DRu,DRv]/DF;


normal=zeros(2,1);
normal(1)=DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);

DR_n=DR*normal;
DR_n=DR_n';

DR_n_2 = DR_n;


temp=Qw(uspan-pu:uspan,1); % 现在对 R_{i,1}的梯度进行处理.

DRu=temp.*((W*DNu-DWu*Nu)*Nv(1))/W2; 
DRv=temp.*(Nu*(DNv(1)*W-DWv*Nv(1)))/W2; 
DR=[DRu,DRv]/DF;
DR_n=DR*normal;
DR_n=DR_n';
DR_n_1 = DR_n;





J=norm(DC,2);

s=(g2(C(1),C(2)) - DR_n_2*u_h_g2(uspan-pu:uspan) - DR_n_1*u_h_g1(uspan-pu:uspan)   ).^2*J;
end
