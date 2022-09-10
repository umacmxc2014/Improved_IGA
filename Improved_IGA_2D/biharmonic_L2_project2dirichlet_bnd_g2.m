function  [g2_h,bnd_dofs_2nd_layers]= biharmonic_L2_project2dirichlet_bnd_g2(ConPts,weights,knotU,pu,knotV,pv,u_grad,g1_h,nurbs_refine)

 
Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
 

nurbs_refine.n
 

Ubreaks=nurbs_refine.UBreaks;
Vbreaks=nurbs_refine.VBreaks;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;




bottom_dofs = nurbs_refine.bottom_dofs;
top_dofs = nurbs_refine.top_dofs;
left_dofs = nurbs_refine.left_dofs;
right_dofs = nurbs_refine.right_dofs;



bnd_dofs_1st_layers = [bottom_dofs,top_dofs,left_dofs,right_dofs];



bnd_dofs_1st_layers = unique(bnd_dofs_1st_layers);


bottom_dofs_2nd_layer = nurbs_refine.bottom_dofs_2nd_layer;
top_dofs_2nd_layer = nurbs_refine.top_dofs_2nd_layer;
left_dofs_2nd_layer = nurbs_refine.left_dofs_2nd_layer;
right_dofs_2nd_layer = nurbs_refine.right_dofs_2nd_layer;


bnd_dofs_2nd_layers = [bottom_dofs_2nd_layer,top_dofs_2nd_layer,left_dofs_2nd_layer, right_dofs_2nd_layer];
bnd_dofs_2nd_layers = unique(bnd_dofs_2nd_layers);



bottom_dofs_2_layers = nurbs_refine.bottom_dofs_2_layers;
top_dofs_2_layers = nurbs_refine.top_dofs_2_layers;
left_dofs_2_layers = nurbs_refine.left_dofs_2_layers ;
right_dofs_2_layers = nurbs_refine.right_dofs_2_layers ;



%% The bottom boundary

[A_bottom,rhs_bottom] = generate_matrix_rhs_g2_bottom(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine);

%%  The top boundary

[A_top,rhs_top] = generate_matrix_rhs_g2_top(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine);

 %%  The left boundary
 
[A_left,rhs_left] = generate_matrix_rhs_g2_left(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine);

%%   The right boundary

[A_right,rhs_right] = generate_matrix_rhs_g2_right(ConPts,weights,knotU,pu,knotV,pv,u_grad,nurbs_refine);
%%

bnd_dofs_2_layers = [bottom_dofs_2_layers,top_dofs_2_layers,left_dofs_2_layers,right_dofs_2_layers];



bnd_dofs_2_layers = unique(bnd_dofs_2_layers);


n_dofs=nurbs_refine.dof;

M = sparse(n_dofs,n_dofs);
F = zeros(n_dofs,1);


M(bottom_dofs_2_layers,bottom_dofs_2_layers) = M(bottom_dofs_2_layers,bottom_dofs_2_layers) + A_bottom;
F(bottom_dofs_2_layers) = F(bottom_dofs_2_layers) + rhs_bottom;
 
M(top_dofs_2_layers,top_dofs_2_layers) = M(top_dofs_2_layers,top_dofs_2_layers) + A_top;
F(top_dofs_2_layers) = F(top_dofs_2_layers) + rhs_top;

M(left_dofs_2_layers,left_dofs_2_layers) = M(left_dofs_2_layers,left_dofs_2_layers)  + A_left;
F(left_dofs_2_layers) = F(left_dofs_2_layers) + rhs_left;


M(right_dofs_2_layers,right_dofs_2_layers) = M(right_dofs_2_layers,right_dofs_2_layers) + A_right;
F(right_dofs_2_layers) = F(right_dofs_2_layers) + rhs_right;

x_non_zero = zeros(n_dofs,1);

x_non_zero(bnd_dofs_1st_layers) = g1_h(bnd_dofs_1st_layers);

F_modify = F - M*x_non_zero;


M(bnd_dofs_1st_layers,:)=0; M(:,bnd_dofs_1st_layers)=0; 
M(bnd_dofs_1st_layers,bnd_dofs_1st_layers)=speye(length(bnd_dofs_1st_layers));
F_modify(bnd_dofs_1st_layers)=0;






M_bnd = M(bnd_dofs_2_layers,bnd_dofs_2_layers);
F_bnd = F_modify(bnd_dofs_2_layers);


u_d_h_bnd =M_bnd\F_bnd;

u_d_h = zeros(n_dofs,1);

u_d_h(bnd_dofs_2_layers) = u_d_h_bnd;

u_d_h = u_d_h + x_non_zero;


% u_d_h_bnd = M_bnd\F_bnd;

g2_h = zeros(n_dofs,1);

g2_h(bnd_dofs_2nd_layers) = u_d_h(bnd_dofs_2nd_layers);




 
 
 


end
 
 
 
 
 
