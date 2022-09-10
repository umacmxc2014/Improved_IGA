function  [g1_h,bnd_dofs] = L2_project2dirichlet_bnd_g1(ConPts,weights,knotU,pu,knotV,pv,g1,nurbs_refine)

 
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
m=nurbs_refine.m;
n=nurbs_refine.n;

 

Ubreaks=nurbs_refine.UBreaks;
Vbreaks=nurbs_refine.VBreaks;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;


[n_conpts_u,n_conpts_v,DIM]=size(ConPts);
% n_conpts_u　代表 原始的NURBS曲面上的　u　方向上的控制点个数.

bottom_ConPts=zeros(DIM,n_conpts_u);
top_ConPts = zeros(DIM,n_conpts_u);
left_ConPts = zeros(DIM,n_conpts_v);
right_ConPts = zeros(DIM,n_conpts_v);


j = 1;
for d=1:DIM
  for i=1:n_conpts_u
    bottom_ConPts(d,i)=ConPts(i,j,d);  % 描述底部边界的控制点信息
end
end
bottom_weights=weights(:,j)'; % 我们要化成行向量

j = n_conpts_v;
for d=1:DIM
  for i=1:n_conpts_u
    top_ConPts(d,i)=ConPts(i,j,d); % 描述顶部边界的控制点信息
end
end
top_weights=weights(:,j)'; % 我们要化成行向量

i = 1;
for d=1:DIM
  for j=1:n_conpts_v
    left_ConPts(d,j)=ConPts(i,j,d); % 描述左边边界的控制点信息
end
end
left_weights = weights(i,:);  % 现在已经是行向量

i = n_conpts_u;
for d=1:DIM
  for j=1:n_conpts_v
    right_ConPts(d,j)=ConPts(i,j,d);  % 描述右端边界的控制点信息
end
end
right_weights = weights(i,:); % 现在已经是行向量




 [A_bottom,rhs_bottom] = L2_project2dirichlet_bnd(bottom_ConPts,knotU,bottom_weights,pu, Ubar,m, Ubreaks, uNoEs, g1);
 [A_top,rhs_top] = L2_project2dirichlet_bnd(top_ConPts,knotU,top_weights,pu, Ubar,m, Ubreaks, uNoEs, g1);
 [A_left,rhs_left] = L2_project2dirichlet_bnd(left_ConPts,knotV,left_weights,pv, Vbar,n, Vbreaks, vNoEs, g1);
 [A_right,rhs_right] = L2_project2dirichlet_bnd(right_ConPts,knotV,right_weights,pv, Vbar,n, Vbreaks, vNoEs, g1);
 
 
 
bottom_dofs = nurbs_refine.bottom_dofs;
top_dofs = nurbs_refine.top_dofs;
left_dofs = nurbs_refine.left_dofs;
right_dofs = nurbs_refine.right_dofs;

bnd_dofs = [bottom_dofs,top_dofs,left_dofs,right_dofs];
bnd_dofs = unique(bnd_dofs);

n_dofs=nurbs_refine.dof;

M = sparse(n_dofs,n_dofs);
F = zeros(n_dofs,1);


M(bottom_dofs,bottom_dofs) = M(bottom_dofs,bottom_dofs) + A_bottom;
F(bottom_dofs) = F(bottom_dofs) + rhs_bottom;
 
M(top_dofs,top_dofs) = M(top_dofs,top_dofs) + A_top;
F(top_dofs) = F(top_dofs) + rhs_top;

M(left_dofs,left_dofs) = M(left_dofs,left_dofs) + A_left;
F(left_dofs) = F(left_dofs) + rhs_left;

M(right_dofs,right_dofs) = M(right_dofs,right_dofs) + A_right;
F(right_dofs) = F(right_dofs) + rhs_right;

M_bnd = M(bnd_dofs,bnd_dofs);
F_bnd  = F(bnd_dofs);

u_d_h_bnd = M_bnd\F_bnd;

g1_h = zeros(n_dofs,1);

g1_h(bnd_dofs) = u_d_h_bnd;



 
 

end




