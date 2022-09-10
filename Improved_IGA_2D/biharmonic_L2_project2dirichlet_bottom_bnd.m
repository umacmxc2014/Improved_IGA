function  [u_h_g1,u_h_g2] = biharmonic_L2_project2dirichlet_bottom_bnd(ConPts,weights,knotU,pu,knotV,pv,Refinement,g1,g2)
% The dimension of ConPts is 



% g1=@(x,y)sin(pi*x+1)*sin(pi*y+1);% u = g1 on $\partial \Omega$.
% g2=@(x,y) sin(pi*x+1)*sin(pi*y+1); % grad u \cdot n =g2  on $\partial \Omega$.



nurbsInfo=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);

 
Ubar=nurbsInfo.Ubar;
Vbar=nurbsInfo.Vbar;
m=nurbsInfo.m;
 

 

Ubreaks=nurbsInfo.UBreaks;
 

uNoEs=nurbsInfo.uNoEs;



[n_conpts_u,~,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

bottom_ConPts=zeros(DIM,n_conpts_u);


j = 1;

for d=1:DIM
  for i=1:n_conpts_u
    bottom_ConPts(d,i)=ConPts(i,j,d);
end
end



bottom_weights=weights(:,j)';




 [A,rhs] = L2_project2dirichlet_bnd(bottom_ConPts,knotU,bottom_weights,pu, Ubar,m, Ubreaks, uNoEs, g1);
 
u_h_g1 = A\rhs;
 ele_dof = pu+1;

coordinate=zeros(uNoEs,2);
element=zeros(uNoEs,ele_dof);
knotSpanIndex=zeros(uNoEs,1);

for i=1:uNoEs
coordinate(i,:)=[Ubreaks(i),Ubreaks(i+1)];
span=findspan(Ubar,pu,Ubreaks(i));
knotSpanIndex(i)=span;
element(i,:)=span-pu:span;
end

n_dofs= m;
A=sparse(2*n_dofs,2*n_dofs);
rhs=zeros(2*n_dofs,1);


np=pu + 1 ;% 每一条边上的数值积分点的个数.
 
 for i=1:uNoEs
	column_idx=element(i,:);
    row_idx = [column_idx, column_idx+m]; 
    ue=coordinate(i,:);
    uspan=knotSpanIndex(i);
    m_row = 2*ele_dof ; 
    n_column = 2*ele_dof + 1 ;
    x1={np,ue(1),ue(2),m_row,n_column};
    s=Gauss_1d(@(u)quad_Ae_Fe_biharmonic(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, uspan, g2, u_h_g1),x1{:});
     
    % size(s)
    A(row_idx,row_idx)=A(row_idx,row_idx )+s(:,1:2*ele_dof);
    % rhs(row_idx)
    rhs(row_idx)=rhs(row_idx)+s(:,end);
 end
 

 

    
    A_21 = A(m+1:end,1:m);  
    
    x_non_zero = 0*rhs;
    
    x_non_zero(1:m) = u_h_g1;
    
    rhs_new = rhs - A*x_non_zero;
    
    
    
    rhs(1:m) = u_h_g1;
    rhs(m+1:end) = rhs(m+1:end) - A_21*u_h_g1;

    
   A(1:m,1:m) = speye(m,m);
   A(1:m,m+1:end) = 0;
   A(m+1:end,1:m) = 0;
   
   rhs_new(1:m) = 0;
    
    x_0 = A\rhs_new;
    
 
    
    x_new = x_0 + x_non_zero;
 % det(A)

 u_h_g2=A\rhs;
 
 


 % u_h_g1 = u_h_g2(1:m);
 
 u_h_g2=u_h_g2(m+1:end);

 
 
  err_g2=0;% err_g2(1)存储的是g2-梯度\Pi \tilde{g}的L2误差。
  % err_g2(２)存储的是g2-梯度\Pi \tilde{g}与的\Pi \tilde{g}的L2内积。
  

  
 for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_biharmonic(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, g2),x1{:});
   err_g2=err_g2+s;
   
 end
 
 
  err_g2 = sqrt(abs(err_g2));
 
 disp('The old one for boundary error is ')
 disp(err_g2)



end
 
 
 


 
 function s=quad_Ae_Fe_biharmonic(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, uspan, g2,u_h_g1)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。

v = 0; % 对于底边， v=0

[~,DF,~,~,~]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);




[n_conpts_u,~,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

bottom_ConPts=zeros(DIM,n_conpts_u);

j = 1;

for d=1:DIM
  for i=1:n_conpts_u
    bottom_ConPts(d,i)=ConPts(i,j,d);
end
end



bottom_weights=weights(:,j)';





[~,~,C,DC]=NurbsCurve(bottom_ConPts,knotU,bottom_weights,pu,u);


normal=zeros(2,1);
normal(1)= DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);

J = norm(DC,2);

% disp('================')

% disp('================')
% if(J<=1.0e-5|| abs(det(DF)) <= 1.0e-5)
%    disp('***********************************')
%    disp('There is something wrong???')
%    disp('***********************************')
% end


Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pv,v,1); %  这时候只需要 v =0 处的基函数及其导数值.

Nu=Uders(1,:)'; DNu=Uders(2,:)';
Nv=Vders(1,:);  DNv=Vders(2,:);

B = Nu; % 当前边上沿着 u 方向上的非0 B-样条基函数组成的列向量;

DB_i1u = DNu*Nv(1); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1v = Nu*DNv(1); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1  = [DB_i1u,DB_i1v]/DF;
DB_i1_n = DB_i1*normal;
D_u_g1_n= (u_h_g1(uspan-pu:uspan)' * DB_i1_n);% g1 的 L2 投影的梯度点乘 n.

DB_i2u = DNu*Nv(2); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2v = Nu*DNv(2); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2  = [DB_i2u,DB_i2v]/DF;
DB_i2_n = DB_i2*normal;
DB_i2_n = DB_i2_n';  % 这个时候是行向量了.

grad_B = [DB_i1u,DB_i1v;DB_i2u,DB_i2v]/DF;
grad_B_n = grad_B*normal;


Ae = grad_B_n*grad_B_n'*J;



Fe=g2(C(1),C(2)) *grad_B_n*J;



s=[Ae,Fe];

 end

 
 
 function s=quad_err_L2_biharmonic(u,Ubar,Vbar, ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, g2)

v = 0;

[~,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

[n_conpts_u,~,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

U_ConPts=zeros(DIM,n_conpts_u);

j = 1;

for d=1:DIM
  for i=1:n_conpts_u
    U_ConPts(d,i)=ConPts(i,j,d);
end
end



U_weights=weights(:,j)';





[~,~,C,DC]=NurbsCurve(U_ConPts,knotU,U_weights,pu,u);

J=norm(DC,2);
normal=zeros(2,1);
normal(1)=DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);



Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,v,1);

Nu=Uders(1,:)'; DNu=Uders(2,:)';
Nv=Vders(1,:); DNv=Vders(2,:);

DB_i1u = DNu*Nv(1); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1v = Nu*DNv(1); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1  = [DB_i1u,DB_i1v]/DF;
DB_i1_n = DB_i1*normal; % 列向量
D_u_g1_n= u_h_g1(uspan-pu:uspan)' * DB_i1_n ;% g1 的 L2 投影的梯度点乘 n.

DB_i2u = DNu*Nv(2); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2v = Nu*DNv(2); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2  = [DB_i2u,DB_i2v]/DF;
DB_i2_n = DB_i2*normal; % 列向量
D_u_g2_n= u_h_g2(uspan-pu:uspan)' * DB_i2_n ;% g1 的 L2 投影的梯度点乘 n.

 s=(g2(C(1),C(2)) - D_u_g1_n - D_u_g2_n  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
% This is used to calculate $\int_{0}^1 (g_2 - \nabla \tilde{g}) \cdot \nabla \tilde{g} d\xi$
%  s=(g2(C(1),C(2)) - D_u_g1_n - D_u_g2_n  ).*(D_u_g1_n + D_u_g2_n)*J;
 

end