function  [A,rhs] = generate_matrix_rhs_g2_left(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine)

% Note that $g_2 = \Delta u \cdot \bm{n}$, on $\partial \Omega$


Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;

n=nurbs_refine.n;
 

 

 
Vbreaks=nurbs_refine.VBreaks;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;


ele_n_dofs = pv+1;

coordinate=zeros(vNoEs,2);
element=zeros(vNoEs,2*ele_n_dofs);
knotSpanIndex=zeros(vNoEs,1);

for j=1:vNoEs
coordinate(j,:)=[Vbreaks(j),Vbreaks(j+1)];
span=findspan(Vbar,pv,Vbreaks(j));
knotSpanIndex(j)=span;
local_idx = 1;
for k=(span-pv):span
    for i1=1:2
element(j,local_idx)=i1+(k-1)*2;
local_idx = local_idx + 1;
    end
end
end

n_dofs = 2*n;
A=sparse(n_dofs,n_dofs);
rhs=zeros(n_dofs,1);


np=pv + 1 ;% 每一条边上的数值积分点的个数.
 
 for j=1:vNoEs
	row=element(j,:);
     
    ve=coordinate(j,:);
   
    m_row = 2*ele_n_dofs ; 
    n_column = 2*ele_n_dofs + 1 ;
    x1={np,ve(1),ve(2),m_row,n_column};
    s=Gauss_1d(@(v)quad_Ae_Fe_biharmonic_left(v,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv,  u_grad),x1{:});
     
    A(row,row)=A(row,row)+s(:,1:2*ele_n_dofs);
    rhs(row)=rhs(row)+s(:,end);
end
 




end
 
 
 


 
 function s=quad_Ae_Fe_biharmonic_left(v,Ubar,Vbar,   ConPts,weights,knotU,knotV, pu,pv,u_grad)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。

u = 0;

[~,DF,~,~,~]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);




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


normal=zeros(2,1);
normal(1)= - DC(2); normal(2)=  DC(1); % 对于 left  boundary, 控制点的顺序为顺时针.
normal=normal/norm(normal);   


g2 = u_grad(C(1),C(2))*normal; % g_2 = grad u \cdot n;




J=norm(DC,2);

Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pv,v,1); %  这时候只需要 v =1 处的基函数及其导数值.


Nv=Vders(1,:);       DNv=Vders(2,:);        %  这些是行向量
Nu=Uders(1,1:2)';  DNu=Uders(2,1:2)';   % 化成列向量


D_B_u = DNu*Nv;  D_B_u = D_B_u(:);
D_B_v = Nu*DNv;  D_B_v = D_B_v(:);

grad_B = [D_B_u,D_B_v]/DF;

grad_B_n = grad_B*normal;


Ae = grad_B_n*grad_B_n'*J;

Fe=g2*grad_B_n*J;



s=[Ae,Fe];

 end