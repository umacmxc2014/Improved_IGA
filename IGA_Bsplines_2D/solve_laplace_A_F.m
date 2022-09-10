function [A,rhs] = solve_laplace_A_F(nurbs_original,nurbs_refine,f)

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


Element=nurbs_refine.Element;
knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;
 
UBreaks=nurbs_refine.UBreaks;   % u 方向上节点向量中的断点.
VBreaks=nurbs_refine.VBreaks;    % v 方向上节点向量中的断点.
 

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;

uNoEs=nurbs_refine.uNoEs;
vNoEs=nurbs_refine.vNoEs;
 

pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;
 

u_np = pu  ;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv ;   
 
[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].

u_ele_basis_funcs_o =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad_o  =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);

u_knot_span_o = zeros(uNoEs,1);
u_jacobian    = zeros(uNoEs, u_np);



basis = zeros(pu_o+1,u_np);
basis_grad = basis;
for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span_o(i) = findspan(knotU_o,pu_o,ue(1));
    for j=1:u_np
     u = uJ*gp_u(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(knotU_o,pu_o,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
     u_jacobian(i,j) = uJ*gw_u(j);
    end
u_ele_basis_funcs_o{i}=basis;
u_ele_basis_grad_o{i}=basis_grad;
end

u_ele_basis_funcs =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_ele_basis_grad  =  cell(uNoEs,1); % zeros(uNoEs,pu+1,u_np);
u_knot_span = zeros(uNoEs,1);
basis = zeros(pu+1,u_np);
basis_grad = basis;
for i=1:uNoEs
	ue = UBreaks(i:i+1);
    uJ=(ue(2)-ue(1))/2;
    u_knot_span(i) = findspan(knotU,pu,ue(1));
    for j=1:u_np
     u = uJ*gp_u(j)+(ue(1)+ue(2))/2;
     Uders=bspbasisDers(knotU,pu,u,1); 
     basis(:,j) = Uders(1,:);
     basis_grad(:,j) = Uders(2,:);
    end
u_ele_basis_funcs{i}=basis;
u_ele_basis_grad{i}=basis_grad;
end



[gp_v,gw_v]=grule(v_np); % The quadrature points and quadrature weights in the [-1,1].

v_ele_basis_funcs_o =  cell(vNoEs,1); % zeros(vNoEs,pv+1,v_np);
v_ele_basis_grad_o  =  cell(vNoEs,1); % zeros(vNoEs,pv+1,v_np);
v_knot_span_o = zeros(vNoEs,1);
v_jacobian    = zeros(vNoEs, v_np);

basis = zeros(pv_o+1,v_np);
basis_grad = basis;
for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    v_knot_span_o(i) = findspan(knotV_o,pv_o,ve(1));
    for j=1:v_np
     v = vJ*gp_v(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(knotV_o,pv_o,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
     v_jacobian(i,j) = vJ*gw_v(j);
    end
v_ele_basis_funcs_o{i}=basis;
v_ele_basis_grad_o{i}=basis_grad;
end

basis = zeros(pv+1,v_np);
basis_grad = basis;
v_ele_basis_funcs =  cell(vNoEs,1); % zeros(uNoEs,pu+1,u_np);
v_ele_basis_grad  =  cell(vNoEs,1); % zeros(uNoEs,pu+1,u_np);
v_knot_span = zeros(vNoEs,1);
for i=1:vNoEs
	ve = VBreaks(i:i+1);
    vJ=(ve(2)-ve(1))/2;
    v_knot_span(i) = findspan(knotV,pv,ve(1));
    for j=1:v_np
     v = vJ*gp_v(j)+(ve(1)+ve(2))/2;
     Vders=bspbasisDers(knotV,pv,v,1); 
     basis(:,j) = Vders(1,:);
     basis_grad(:,j) = Vders(2,:);
    end
v_ele_basis_funcs{i}=basis;
v_ele_basis_grad{i}=basis_grad;
end








n_gps = u_np*v_np;


n_ele_dofs = (pu+1)*(pv+1);

Ae =zeros(n_ele_dofs,n_ele_dofs);
Fe = zeros(n_ele_dofs,1);

% A = sparse(n_dofs,n_dofs);
rhs = zeros(n_dofs,1);

DF = cell(n_gps,1);
Jacobian = zeros(n_gps,1);

value = zeros(NoEs*n_ele_dofs*n_ele_dofs,1);
row_index = value;
column_index = value;
global_index = 1;


    for j=1:vNoEs
        for i=1:uNoEs
            Ae = 0*Ae; Fe =0*Fe; 
            e = i + (j-1)*uNoEs ;
            row = Element(e,:);

            u_basis_i_o   =  u_ele_basis_funcs_o{i};
            v_basis_j_o   =  v_ele_basis_funcs_o{j};
             
            uv_basis_ij_o = kron(v_basis_j_o,u_basis_i_o);
            
            
            Du_basis_i_o =  u_ele_basis_grad_o{i};
            Du_v_basis_ij_o = kron(v_basis_j_o,Du_basis_i_o);
             
            
            Dv_basis_j_o =  v_ele_basis_grad_o{j};
            Dv_u_basis_ij_o = kron(Dv_basis_j_o,u_basis_i_o);
             
           
            
            
            uspan_o  =  u_knot_span_o(i);
            vspan_o  =  v_knot_span_o(j);
            
            
            w = weights_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o);
            w = w(:);
            w_funcs_o = w'*uv_basis_ij_o; % 现在算出了权重函数在 u_np*v_np*w_np 个积分点处的值
            
            
            P_x = ConPts_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,1);
            P_x = P_x(:);
            w_P_x = (w.*P_x)';
            F_x = (w_P_x* uv_basis_ij_o)./w_funcs_o;
            
            
            P_y = ConPts_o(uspan_o-pu_o:uspan_o,vspan_o-pv_o:vspan_o,2);
            P_y = P_y(:);
            w_P_y = (w.*P_y)';
            F_y = (w_P_y* uv_basis_ij_o)./w_funcs_o;
            

            
            F=[F_x;F_y]; % 2*n_gps
            
            
            Dw_u_o   = w'*Du_v_basis_ij_o; 
            Dw_v_o   = w'*Dv_u_basis_ij_o; 
            
            
        
            DF_x_u  = (w_P_x*Du_v_basis_ij_o - F_x.*Dw_u_o )./w_funcs_o;
            DF_x_v  = (w_P_x*Dv_u_basis_ij_o - F_x.*Dw_v_o )./w_funcs_o;
    

            DF_y_u  = (w_P_y*Du_v_basis_ij_o - F_y.*Dw_u_o )./w_funcs_o;
            DF_y_v  = (w_P_y*Dv_u_basis_ij_o - F_y.*Dw_v_o )./w_funcs_o;
             
            

            
            
            u_gw = u_jacobian(i,:);         u_gw = u_gw(:);
            v_gw = v_jacobian(j,:);         v_gw  = v_gw(:);
             
            uv_gw = kron(v_gw,u_gw); 
           
            
            for k1=1:n_gps
            DF{k1}=[DF_x_u(k1)  DF_x_v(k1);... 
                            DF_y_u(k1)  DF_y_v(k1)];
            Jacobian(k1) = abs(det(DF{k1}))*uv_gw(k1);
            end
            

            
            
            
            
            u_basis_i   =  u_ele_basis_funcs{i};
            v_basis_j   =  v_ele_basis_funcs{j};
             
            uv_basis_ij = kron(v_basis_j,u_basis_i);
            
            
            Du_basis_i  = u_ele_basis_grad{i};
            Du_v_basis_ij = kron(v_basis_j,Du_basis_i);
           
            
            Dv_basis_j  = v_ele_basis_grad{j};
            Dv_u_basis_ij = kron(Dv_basis_j,u_basis_i);
           
            
           
           
            
            for k1=1:n_gps
                basis_grad = [Du_v_basis_ij(:,k1),Dv_u_basis_ij(:,k1)]/DF{k1};
                Ae = Ae + basis_grad*basis_grad'*Jacobian(k1);
                Fe = Fe + f(F(1,k1),F(2,k1))*uv_basis_ij(:,k1)*Jacobian(k1);
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