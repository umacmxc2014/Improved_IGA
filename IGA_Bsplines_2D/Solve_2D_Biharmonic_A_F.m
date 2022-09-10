function [A,rhs] = Solve_2D_Biharmonic_A_F(nurbs_original,nurbs_refine,f) 


ConPts_o = nurbs_original.ConPts;
weights_o = nurbs_original.weights;
knotU = nurbs_original.knotU;
knotV = nurbs_original.knotV;
pu_o   =  nurbs_original.pu;
pv_o   =  nurbs_original.pv;

DIM = 2;


Element=nurbs_refine.Element;
Coordinate=nurbs_refine.Coordinate;

Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
pu=nurbs_refine.pu;
pv=nurbs_refine.pv;
 

 

NoEs=nurbs_refine.NoEs;
n_dofs=nurbs_refine.n_dofs;


 

u_np = pu + 1 ;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv + 1;   
 
[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].

gp_v = gp_u;
gw_v = gw_u;

ele_n_dofs = (pu+1)*(pv+1);

rhs = zeros(n_dofs,1);

Ae = zeros(ele_n_dofs,ele_n_dofs);

Fe = zeros(ele_n_dofs,1);


value = zeros(NoEs*ele_n_dofs*ele_n_dofs,1);
row_index = value;
column_index = value;
global_index = 1;


I_2 = eye(DIM);% 2 by 2 identity matrix.

D_DF_x = zeros(DIM,DIM);



for e=1:NoEs % Loop for each element;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);

     
    
    Ae = 0*Ae; Fe = 0*Fe;
    
    for j=1:v_np
            vJ = (ve(2) - ve(1))/2;
            v = vJ*gp_v(j) + (ve(1)+ve(2))/2;
            Vders=bspbasisDers(Vbar,pv,v,2);
            Nv=Vders(1,:); DNv=Vders(2,:);  D2Nv=Vders(3,:);
        for i=1:u_np
            uJ = (ue(2) - ue(1))/2;
            u = uJ*gp_u(i)+(ue(1)+ue(2))/2;
            Uders=bspbasisDers(Ubar,pu,u,2);
            Nu=Uders(1,:)';DNu=Uders(2,:)';  D2Nu=Uders(3,:)';
            [~,~, ~,F,DF,D2F]=NurbsSurfaceDers(ConPts_o,knotU,knotV,weights_o,pu_o,u,pv_o,v);
            
            inv_DF=DF\I_2;

            
            D_DF_x(1,1) = [D2F(1,1,1), D2F(1,2,1)] *inv_DF(:,1);
            D_DF_x(1,2) = [D2F(1,1,2), D2F(1,2,2)]*inv_DF(:,1);
            D_DF_x(2,1) = [D2F(2,1,1), D2F(2,2,1)]*inv_DF(:,1);
            D_DF_x(2,2) = [D2F(2,1,2), D2F(2,2,2)]*inv_DF(:,1);

            D_inv_DF_x =  - inv_DF*D_DF_x*inv_DF;


            D_DF_y = zeros(DIM,DIM);
            D_DF_y(1,1) = [D2F(1,1,1), D2F(1,2,1)] *inv_DF(:,2);
            D_DF_y(1,2) = [D2F(1,1,2), D2F(1,2,2)]*inv_DF(:,2);
            D_DF_y(2,1) = [D2F(2,1,1), D2F(2,2,1)]*inv_DF(:,2);
            D_DF_y(2,2) = [D2F(2,1,2), D2F(2,2,2)]*inv_DF(:,2);

            D_inv_DF_y = - inv_DF*D_DF_y*inv_DF;

            
            J=abs(det(DF))*gw_u(i)*uJ*gw_v(j)*vJ;
           
            B = Nu*Nv;  B=B(:);
            

            Fe =Fe + f(F(1),F(2))*B*J;

                  DBu =DNu*Nv; DBu =DBu(:);      
                  DBv =Nu*DNv; DBv =DBv(:);     
                  

                  D2_B_uu  = D2Nu * Nv ; D2_B_uu = D2_B_uu(:);      
                  D2_B_vv  = Nu * D2Nv ; D2_B_vv = D2_B_vv(:);       
                  D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:);     

                 
                 D2_B = [ [D2_B_uu,D2_B_uv];[D2_B_uv,D2_B_vv ] ];
                 
                 % D2_B存储的是$\grad \hat{u}$的三个分量的梯度,合在一块就是三乘三的矩阵.
                 
                 DB_x = D2_B*inv_DF(:,1);  % DB_x代表的是 $\grad \hat{u}$关于x_1的偏导数
                 DB_x = reshape(DB_x,ele_n_dofs,DIM);
                 
                 DB_y = D2_B*inv_DF(:,2); % DB_x代表的是 $\grad \hat{u}$关于x_2的偏导数
                 DB_y = reshape(DB_y,ele_n_dofs,DIM);
                 

                 
                 
               D2_B_xx=DB_x*inv_DF(:,1) +  [DBu,DBv]*D_inv_DF_x(:,1);
               D2_B_yy=DB_y*inv_DF(:,2) +  [DBu,DBv]*D_inv_DF_y(:,2);
           


                 
               Laplacian_B = D2_B_xx + D2_B_yy  ;
                                       



                Ae = Ae + Laplacian_B*Laplacian_B'*J; % This is the element stiffness matrix for the  biharmonic equation.

        end
    end
    
  
            
            for i1=1:ele_n_dofs
               for j1=1:ele_n_dofs
                    row_index(global_index) = row(i1);
                    column_index(global_index) = row(j1);
                    value(global_index) = value(global_index) + Ae(i1,j1);
                   global_index = global_index + 1;
               end
            end
           
            
            rhs(row) = rhs(row)+Fe;
            
end



A =sparse(row_index,column_index,value,n_dofs,n_dofs); 




end
