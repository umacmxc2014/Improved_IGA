function [A,rhs] = Solve_Biharmonic_A_Mass(nurbs_original,nurbs_refine,f) 


DIM = 2;

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


Element=nurbs_refine.Element;
knotU=nurbs_refine.Ubar;
 
 
UBreaks=nurbs_refine.UBreaks;   % u 方向上节点向量中的断点.
 
 

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;
 



pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
 

u_np = pu +1 ;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv +1;   
 
[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].



u_ele_basis_funcs =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_der  =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_2ders = cell(uNoEs,1);

u_knot_span = zeros(uNoEs,1);
u_quad_pnts = zeros(uNoEs,u_np);

basis = zeros(pu+1,u_np);
basis_grad = basis;
basis_hessian = basis;
u_jocabian = zeros(uNoEs,u_np);

for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span(i) = findspan(knotU,pu,ue(1));
    for j=1:u_np
     u = uJ*gp_u(j)+(ue(1)+ue(2))/2;
     u_quad_pnts(i,j) = u;
     Uders=bspbasisDers(knotU,pu,u,2); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
     basis_hessian(:,j) = Uders(3,:);
     u_jocabian(i,j) = uJ*gw_u(j);
    end
u_ele_basis_funcs{i}=basis;
u_ele_basis_der{i}=basis_grad;
u_ele_basis_2ders{i} = basis_hessian;
end


v_ele_basis_funcs = u_ele_basis_funcs;
v_ele_basis_der  = u_ele_basis_der;
v_quad_pnts = u_quad_pnts;
v_ele_basis_2ders = u_ele_basis_2ders;
v_jocabian = u_jocabian;



n_ele_dofs = (pu+1)*(pv+1);

Ae =zeros(n_ele_dofs,n_ele_dofs);
Fe = zeros(n_ele_dofs,1);

Me = zeros(n_ele_dofs,n_ele_dofs);



% A = sparse(n_dofs,n_dofs);
rhs = zeros(n_dofs,1);

value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
row_index = value;
column_index = value;
global_index = 1;


M_value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);


 I_2 = eye(2);% 2 by 2 identity matrix.
 
 D_DF_x = zeros(DIM,DIM);
 
 D_DF_y = zeros(DIM,DIM);

    for j=1:vNoEs
        for i=1:uNoEs
            Ae = 0*Ae; Fe =0*Fe; 
            e = i + (j-1)*uNoEs ;
            row = Element(e,:);
            
            Nu_i=u_ele_basis_funcs{i}; DNu_i=u_ele_basis_der{i}; D2Nu_i=u_ele_basis_2ders{i};
            Nv_j=v_ele_basis_funcs{j}; DNv_j =v_ele_basis_der{j}; D2Nv_j=v_ele_basis_2ders{j};
            
            for i1=1:u_np
                for j1=1:v_np
                    v = v_quad_pnts(j,j1);
                    u = u_quad_pnts(i,i1);
                    Nu=Nu_i(:,i1); DNu=DNu_i(:,i1);    D2Nu=D2Nu_i(:,i1);
                    Nv=Nv_j(:,j1)'; DNv=DNv_j(:,j1)';   D2Nv=D2Nv_j(:,j1)';
                    
                    [~,~, ~,F,DF,D2F]=NurbsSurfaceDers(ConPts_o,knotU_o,knotV_o,weights_o,pu_o,u,pv_o,v);
                    
                   

                    inv_DF=DF\I_2;

                    
                    D_DF_x(1,1) = [D2F(1,1,1), D2F(1,2,1)] *inv_DF(:,1);
                    D_DF_x(1,2) = [D2F(1,1,2), D2F(1,2,2)]*inv_DF(:,1);
                    D_DF_x(2,1) = [D2F(2,1,1), D2F(2,2,1)]*inv_DF(:,1);
                    D_DF_x(2,2) = [D2F(2,1,2), D2F(2,2,2)]*inv_DF(:,1);

                    D_inv_DF_x =  - inv_DF*D_DF_x*inv_DF;


                   
                   D_DF_y(1,1) = [D2F(1,1,1), D2F(1,2,1)] *inv_DF(:,2);
                   D_DF_y(1,2) = [D2F(1,1,2), D2F(1,2,2)]*inv_DF(:,2);
                   D_DF_y(2,1) = [D2F(2,1,1), D2F(2,2,1)]*inv_DF(:,2);
                   D_DF_y(2,2) = [D2F(2,1,2), D2F(2,2,2)]*inv_DF(:,2);

                   D_inv_DF_y = - inv_DF*D_DF_y*inv_DF;
                
                   J=abs(det(DF))*u_jocabian(i,i1)*v_jocabian(j,j1);
                   
                   
                   
                   B = Nu*Nv;  B=B(:);
                   
                   Fe = Fe + f(F(1),F(2))*B*J;


                  DBu=DNu*Nv; DBu=DBu(:);
                  DBv=Nu*DNv; DBv=DBv(:);

                  D2_B_uu= D2Nu * Nv ; D2_B_uu = D2_B_uu(:);
                  D2_B_vv= Nu * D2Nv ; D2_B_vv = D2_B_vv(:);
                  D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:); 



                Laplacian_B=D2_B_uu*sum(inv_DF(1,:).^2) + D2_B_uv*2*(inv_DF(1,1)*inv_DF(2,1)+inv_DF(1,2)*inv_DF(2,2))+ ...
                                     D2_B_vv*sum(inv_DF(2,:).^2) + ...
                                     DBu*(D_inv_DF_x(1,1) + D_inv_DF_y(1,2)) + DBv*(D_inv_DF_x(2,1) + D_inv_DF_y(2,2));



                Ae = Ae + Laplacian_B*Laplacian_B'*J; % This is the element stiffness matrix for the  biharmonic equation.
 
                end
            end
            
            
       
            

             
            for i1=1:n_ele_dofs
               for j1=1:n_ele_dofs
                    row_index(global_index) = row(i1);
                    column_index(global_index) = row(j1);
                    value(global_index) = value(global_index) + Ae(i1,j1);
                   global_index = global_index + 1;
               end
            end
           
            
            % A(row,row) = A(row,row) + Ae;
            rhs(row) = rhs(row)+Fe;
                                           

        end
    end
    

    A =sparse(row_index,column_index,value,n_dofs,n_dofs);     

end