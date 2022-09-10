function  err_g1= bnd_g1_error(ConPts,weights,knotU,pu,knotV,pv,g1,g1_h,nurbs_refine)


 
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
 

Ubreaks=nurbs_refine.UBreaks;
Vbreaks=nurbs_refine.VBreaks;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;




 

np = pu+2;
ele_n_dofs = pu+1;


err_g1 = 0;% err_g2(1)存储的是g2-梯度\Pi \tilde{g}的L2误差;

err_g1_v_0 = 0;
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


for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_g1_bottom(u,Ubar, ConPts,weights,knotU,pu, uspan, g1_h_bottom, g1),x1{:});
   err_g1=err_g1+s;
   
   err_g1_v_0 = err_g1_v_0 +s;
   
end
 


top_dofs = nurbs_refine.top_dofs;

g1_h_top = g1_h(top_dofs);



err_g1_v_1 = 0;

for i=1:uNoEs
 
   ue=coordinate(i,:);
   uspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ue(1),ue(2),m_row,n_column};
   s=Gauss_1d(@(u)quad_err_L2_g1_top(u,Ubar, ConPts,weights,knotU, pu, uspan, g1_h_top, g1),x1{:});
   err_g1=err_g1+s;
   
   err_g1_v_1 = err_g1_v_1 +s;
   
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


err_g1_u_0 = 0;

for i=1:vNoEs
 
   ve=coordinate(i,:);
   vspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ve(1),ve(2),m_row,n_column};
   s=Gauss_1d(@(v)quad_err_L2_g1_left(v,Vbar,  ConPts,weights,knotV,pv, vspan, g1_h_left, g1),x1{:});
   err_g1=err_g1+s;
   
   err_g1_u_0 = err_g1_u_0 +s;
   
end


right_dofs = nurbs_refine.right_dofs;

g1_h_right = g1_h(right_dofs);



err_g1_u_1 = 0;

for i=1:vNoEs
 
   ve=coordinate(i,:);
   vspan=knotSpanIndex(i);
   m_row = 1;
   n_column = 1;
   x1={np,ve(1),ve(2),m_row,n_column};
   s=Gauss_1d(@(v)quad_err_L2_g1_right(v,Vbar,  ConPts,weights,knotV, pv, vspan, g1_h_right, g1),x1{:});
   err_g1=err_g1+s;
   
   err_g1_u_1 = err_g1_u_1 +s;
   
end
 
 
 err_g1 = sqrt(err_g1);
 
 % err_g2_v_0 = sqrt(err_g2_v_0)
 % err_g2_v_1 = sqrt(err_g2_v_1)
 % err_g2_u_0 = sqrt(err_g2_u_0)
 % err_g2_u_1 = sqrt(err_g2_u_1)
 disp('The new one for boundary error is ')

 disp(err_g1)
 % disp(err_g2_v_0)
 % disp(err_g2_v_1)
 

 


end
 
 
 
 
 
 function s=quad_err_L2_g1_bottom(u,Ubar, ConPts,weights,knotU, pu,uspan, u_h_g1,g1)

 
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



Uders=bspbasisDers(Ubar,pu,u,0);


Nu=Uders(1,:); 


 s=(g1(C(1),C(2)) - Nu*u_h_g1(uspan-pu:uspan)  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end

 
function s=quad_err_L2_g1_top(u,Ubar, ConPts,weights,knotU, pu,uspan, u_h_g1,g1)

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



Uders=bspbasisDers(Ubar,pu,u,0);


Nu=Uders(1,:); 


 s=(g1(C(1),C(2)) - Nu*u_h_g1(uspan-pu:uspan)  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$ 
 
end
 



function s=quad_err_L2_g1_left(v,Vbar, ConPts,weights,knotV, pv, vspan, u_h_g1, g1)


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


Vders=bspbasisDers(Vbar,pv,v,0);

Nv=Vders(1,:);


 s=(g1(C(1),C(2)) - Nv*u_h_g1(vspan-pv:vspan)  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end
 
 
 function s=quad_err_L2_g1_right(v,Vbar, ConPts,weights,knotV, pv, vspan, u_h_g1, g1)


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


Vders=bspbasisDers(Vbar,pv,v,0);

Nv=Vders(1,:);


 s=(g1(C(1),C(2)) - Nv*u_h_g1(vspan-pv:vspan)  ).^2*J;% This is used to calculate $\| g_2 - \nabla \tilde{g} \|_{0,\Gamma_D}$
 
 
 end
 
 
 
 
 