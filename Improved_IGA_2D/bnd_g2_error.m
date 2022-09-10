function  err_g2= bnd_g2_error(ConPts,weights,knotU,pu,knotV,pv,u_grad,g1_h,g2_h,nurbs_refine)


 
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
 

Ubreaks=nurbs_refine.UBreaks;
Vbreaks=nurbs_refine.VBreaks;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;





bottom_dofs_2nd_layer = nurbs_refine.bottom_dofs_2nd_layer;
top_dofs_2nd_layer = nurbs_refine.top_dofs_2nd_layer;
left_dofs_2nd_layer = nurbs_refine.left_dofs_2nd_layer;
right_dofs_2nd_layer = nurbs_refine.right_dofs_2nd_layer;

 

np = pu+2;
ele_n_dofs = pu+1;


err_g2 = 0;% err_g2(1)存储的是g2-梯度\Pi \tilde{g}的L2误差;

err_g2_v_0 = 0;
coordinate=zeros(uNoEs,2);
element=zeros(uNoEs,ele_n_dofs);
knotSpanIndex=zeros(uNoEs,1);

for i=1:uNoEs
coordinate(i,:)=[Ubreaks(i),Ubreaks(i+1)];
span=findspan(Ubar,pu,Ubreaks(i));
knotSpanIndex(i)=span;
element(i,:)=span-pu:span;
end



bottom_dofs = nurbs_refine.bottom_dofs;

g1_h_bottom = g1_h(bottom_dofs);

g2_h_bottom = g2_h(bottom_dofs_2nd_layer);

% g2_h_bottom = A_bottom\rhs_bottom;

for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_biharmonic_bottom(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, uspan, g1_h_bottom, g2_h_bottom, u_grad),x1{:});
   err_g2=err_g2+s;
   
   err_g2_v_0 = err_g2_v_0 +s;
   
end
 


top_dofs = nurbs_refine.top_dofs;

g1_h_top = g1_h(top_dofs);

g2_h_top = g2_h(top_dofs_2nd_layer);

% g2_h_top = A_top\rhs_top;

err_g2_v_1 = 0;

for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_biharmonic_top(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, uspan, g1_h_top, g2_h_top, u_grad),x1{:});
   err_g2=err_g2+s;
   
   err_g2_v_1 = err_g2_v_1 +s;
   
end
 


coordinate=zeros(vNoEs,2);
element=zeros(vNoEs,ele_n_dofs);
knotSpanIndex=zeros(vNoEs,1);

for i=1:vNoEs
coordinate(i,:)=[Vbreaks(i),Vbreaks(i+1)];
span=findspan(Vbar,pv,Vbreaks(i));
knotSpanIndex(i)=span;
element(i,:)=span-pv:span;
end



left_dofs = nurbs_refine.left_dofs;

g1_h_left = g1_h(left_dofs);

g2_h_left = g2_h(left_dofs_2nd_layer);

% g2_h_left = A_left\rhs_left;

err_g2_u_0 = 0;

for i=1:vNoEs
 
   ve=coordinate(i,:);
   vspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ve(1),ve(2),m_row,n_column};
   s=Gauss_1d(@(v)quad_err_L2_biharmonic_left(v,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, vspan, g1_h_left, g2_h_left, u_grad),x1{:});
   err_g2=err_g2+s;
   
   err_g2_u_0 = err_g2_u_0 +s;
   
end


right_dofs = nurbs_refine.right_dofs;

g1_h_right = g1_h(right_dofs);

g2_h_right = g2_h(right_dofs_2nd_layer);

% g2_h_right = A_right\rhs_right;

err_g2_u_1 = 0;

for i=1:vNoEs
 
   ve=coordinate(i,:);
   vspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ve(1),ve(2),m_row,n_column};
   s=Gauss_1d(@(v)quad_err_L2_biharmonic_right(v,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, vspan, g1_h_right, g2_h_right, u_grad),x1{:});
   err_g2=err_g2+s;
   
   err_g2_u_1 = err_g2_u_1 +s;
   
end
 
 
 err_g2 = sqrt(err_g2);
 
 % err_g2_v_0 = sqrt(err_g2_v_0)
 % err_g2_v_1 = sqrt(err_g2_v_1)
 % err_g2_u_0 = sqrt(err_g2_u_0)
 % err_g2_u_1 = sqrt(err_g2_u_1)
 disp('The new one for boundary error is ')

 disp(err_g2)
 % disp(err_g2_v_0)
 % disp(err_g2_v_1)
 

 


end
 
 
 
 
 
 function s=quad_err_L2_biharmonic_bottom(u,Ubar,Vbar, ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, u_grad)

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

g2 = u_grad(C(1),C(2))*normal;

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

 s=(g2 - D_u_g1_n - D_u_g2_n  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end

 
function s=quad_err_L2_biharmonic_top(u,Ubar,Vbar, ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, u_grad)

v = 1;

[~,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

[n_conpts_u,n_conpts_v,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

U_ConPts=zeros(DIM,n_conpts_u);

j = n_conpts_v;

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

g2 = u_grad(C(1),C(2))*normal;

Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,v,1);

Nu=Uders(1,:)'; DNu=Uders(2,:)';
Nv=Vders(1,:); DNv=Vders(2,:);

DB_i1u = DNu*Nv(end); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1v = Nu*DNv(end); % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1  = [DB_i1u,DB_i1v]/DF;
DB_i1_n = DB_i1*normal; % 列向量
D_u_g1_n= u_h_g1(uspan-pu:uspan)' * DB_i1_n ;% g1 的 L2 投影的梯度点乘 n.

DB_i2u = DNu*Nv(end-1); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2v = Nu*DNv(end-1); % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2  = [DB_i2u,DB_i2v]/DF;
DB_i2_n = DB_i2*normal; % 列向量
D_u_g2_n= u_h_g2(uspan-pu:uspan)' * DB_i2_n ;% g1 的 L2 投影的梯度点乘 n.

 s=(g2 - D_u_g1_n - D_u_g2_n  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
end
 


 function s=quad_err_L2_biharmonic_left(v,Ubar,Vbar, ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, u_grad)

u = 0;

[~,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

[~,n_conpts_v,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

V_ConPts=zeros(DIM,n_conpts_v);

i = 1;

for d=1:DIM
  for j=1:n_conpts_v
    V_ConPts(d,j)=ConPts(i,j,d);
end
end



V_weights=weights(i,:);





[~,~,C,DC]=NurbsCurve(V_ConPts,knotV,V_weights,pv,v);

J=norm(DC,2);
normal=zeros(2,1);
normal(1)=DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);

g2 = u_grad(C(1),C(2))*normal;

Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,v,1);

Nu=Uders(1,:); DNu=Uders(2,:);
Nv=Vders(1,:)'; DNv=Vders(2,:)';

DB_i1u = DNu(1)*Nv; % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1v = Nu(1)*DNv; % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1  = [DB_i1u,DB_i1v]/DF;
DB_i1_n = DB_i1*normal; % 列向量
D_u_g1_n= u_h_g1(uspan-pu:uspan)' * DB_i1_n ;% g1 的 L2 投影的梯度点乘 n.

DB_i2u = DNu(2)*Nv; % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2v = Nu(2)*DNv; % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2  = [DB_i2u,DB_i2v]/DF;
DB_i2_n = DB_i2*normal; % 列向量
D_u_g2_n= u_h_g2(uspan-pu:uspan)' * DB_i2_n ;% g1 的 L2 投影的梯度点乘 n.

 s=(g2 - D_u_g1_n - D_u_g2_n  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end
 
 
 function s=quad_err_L2_biharmonic_right(v,Ubar,Vbar, ConPts,weights,knotU,knotV, pu,pv, uspan, u_h_g1, u_h_g2, u_grad)

u = 1;

[~,DF]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);

[n_conpts_u,n_conpts_v,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

V_ConPts=zeros(DIM,n_conpts_v);

i = n_conpts_u;

for d=1:DIM
  for j=1:n_conpts_v
    V_ConPts(d,j)=ConPts(i,j,d);
end
end



V_weights=weights(i,:);





[~,~,C,DC]=NurbsCurve(V_ConPts,knotV,V_weights,pv,v);

J=norm(DC,2);
normal=zeros(2,1);
normal(1)=DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);

g2 = u_grad(C(1),C(2))*normal;

Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pu,v,1);

Nu=Uders(1,:); DNu=Uders(2,:);
Nv=Vders(1,:)'; DNv=Vders(2,:)';

DB_i1u = DNu(end)*Nv; % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1v = Nu(end)*DNv; % 当前边上的全体非0的 N_{i,p}(u)N_{1,q}(v)关于u的偏导数
DB_i1  = [DB_i1u,DB_i1v]/DF;
DB_i1_n = DB_i1*normal; % 列向量
D_u_g1_n= u_h_g1(uspan-pu:uspan)' * DB_i1_n ;% g1 的 L2 投影的梯度点乘 n.

DB_i2u = DNu(end-1)*Nv; % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2v = Nu(end-1)*DNv; % 当前边上的全体非0的 N_{i,p}(u)N_{2,q}(v)关于u的偏导数
DB_i2  = [DB_i2u,DB_i2v]/DF;
DB_i2_n = DB_i2*normal; % 列向量
D_u_g2_n= u_h_g2(uspan-pu:uspan)' * DB_i2_n ;% g1 的 L2 投影的梯度点乘 n.

 s=(g2 - D_u_g1_n - D_u_g2_n  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end
 
 