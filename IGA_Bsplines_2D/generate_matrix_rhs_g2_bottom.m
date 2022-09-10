function  [A,rhs] = generate_matrix_rhs_g2_bottom(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine)

% Note that $g_2 = \Delta u \cdot \bm{n}$, on $\partial \Omega$

 
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
m=nurbs_refine.m;
 

 

Ubreaks=nurbs_refine.UBreaks;
 

uNoEs=nurbs_refine.uNoEs;



 ele_n_dofs = pu+1; % 一维单元上非0的基函数个数

coordinate=zeros(uNoEs,2);
element=zeros(uNoEs,2*ele_n_dofs);
knotSpanIndex=zeros(uNoEs,1);

for i=1:uNoEs
coordinate(i,:)=[Ubreaks(i),Ubreaks(i+1)];
span=findspan(Ubar,pu,Ubreaks(i));
knotSpanIndex(i)=span;
local_idx = 1;
for j1=1:2
    for i1=(span-pu):span
element(i,local_idx)=i1+(j1-1)*m;
local_idx = local_idx +1;
    end
end
end




n_dofs= 2*m;
A=sparse(n_dofs,n_dofs);
rhs=zeros(n_dofs,1);


np=pu + 1 ;% 每一条边上的数值积分点的个数.
 
 for i=1:uNoEs
	row_idx=element(i,:);
   
     
    ue=coordinate(i,:);
    
    m_row = 2*ele_n_dofs ; 
    n_column = 2*ele_n_dofs + 1 ;
    x1={np,ue(1),ue(2),m_row,n_column};
    s=Gauss_1d(@(u)quad_Ae_Fe_biharmonic(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv, u_grad),x1{:});
     
    A(row_idx,row_idx)=A(row_idx,row_idx)+s(:,1:2*ele_n_dofs);
    rhs(row_idx)=rhs(row_idx)+s(:,end);
 end
 

 
   

end
 
 
 


 
 function s=quad_Ae_Fe_biharmonic(u,Ubar,Vbar,  ConPts,weights,knotU,knotV, pu,pv,u_grad)

% 把非齐次狄利克雷边界上的边界函数 L^2 投影到边界上的NURBS空间。
%　这里计算L^2投影时，每一条边上所需要的质量矩阵和右端项。

v = 0;

[F,DF,~,~,~]=NurbsSurface(ConPts,weights,knotU,pu,u,knotV,pv,v);



disp('The point is F')
disp(F)

disp('DF is')

disp(DF)

disp('J is')
disp(det(DF))

I_2 = eye(2);

nn=det(DF)*(DF\I_2)'*[0;-1];

nn/norm(nn)


 


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


normal=zeros(2,1); % The unit outward normal vector at point $(x,y) = C(\xi)$.
normal(1)= DC(2); normal(2)= - DC(1);
normal=normal/norm(normal);


disp('DC is')

disp(DC)

disp('n is ')
disp(normal)

g2 = u_grad(C(1),C(2))*normal; % g_2 = grad u \cdot n;

J=norm(DC,2);




Uders=bspbasisDers(Ubar,pu,u,1);
Vders=bspbasisDers(Vbar,pv,v,1); %  这时候只需要 v =0 处的基函数及其导数值.

Nu=Uders(1,:)';     DNu=Uders(2,:)';
Nv=Vders(1,1:2);  DNv=Vders(2,1:2);

D_B_u = DNu*Nv;  D_B_u = D_B_u(:);
D_B_v = Nu*DNv;  D_B_v = D_B_v(:);

grad_B = [D_B_u,D_B_v]/DF;

grad_B_n = grad_B*normal;


Ae = grad_B_n*grad_B_n'*J;

Fe=g2*grad_B_n*J;



s=[Ae,Fe];

 end

 
