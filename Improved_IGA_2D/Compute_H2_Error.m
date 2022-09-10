function err=Compute_H2_Error(nurbs_original,nurbs_refine,Uh,u_Exact, u_Grad,Laplacian_u,D2_u_xy)

% This script is used to compute the $L^2$ norm, $H^1$ semi-norm and $H^2$
% semi-norm  error between $u$ and the IGA solution $u_h$.



DIM = 2;

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;




Element=nurbs_refine.Element;
Coordinate = nurbs_refine.Coordinate;


Ubar=nurbs_refine.Ubar;
Vbar=nurbs_refine.Vbar;
 

 

NoEs=nurbs_refine.NoEs;





 

u_np = pu + 1;   % It seems that if pu >=3, u_np = pu is enough!!!
v_np = pv + 1;   
 
[gp_u,gw_u]=grule(u_np); % The quadrature points and quadrature weights in the [-1,1].

gp_v = gp_u;
gw_v = gw_u;



L2_err = 0;  
H1_err = 0;
H2_err = 0;



I_2 = eye(DIM);% 2 by 2 identity matrix.



D_DF_x = zeros(DIM,DIM);

D_DF_y = zeros(DIM,DIM);


for e=1:NoEs % Loop for each element;
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);

    U_ij = Uh(row);
    U_ij = U_ij';

    
    for j=1:v_np
            vJ = (ve(2) - ve(1))/2;
            v = vJ*gp_v(j) + (ve(1)+ve(2))/2;
            Vders=bspbasisDers(Vbar,pv,v,2);
            Nv=Vders(1,:); DNv=Vders(2,:); D2Nv =Vders(3,:);
        for i=1:u_np
            uJ = (ue(2) - ue(1))/2;
            u = uJ*gp_u(i)+(ue(1)+ue(2))/2;
            Uders=bspbasisDers(Ubar,pu,u,2);
            Nu=Uders(1,:)';DNu=Uders(2,:)'; D2Nu=Uders(3,:)';
            
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
            
            J=abs(det(DF))*gw_u(i)*uJ*gw_v(j)*vJ;
           
            B= Nu*Nv;
            B= B(:);
            
            if (ue(1)==0.5 && ve(2)==100)
            disp('***********************************')
            disp('The point is ')
            disp([F(1),F(2)])
            disp('ue and ve are')
            disp(ue)
            disp(ve)
            disp('The error is ')
            disp(u_Exact(F(1),F(2)) - U_ij*B)
            
            disp('***********************************')
            end
            
            L2_err = L2_err + (u_Exact(F(1),F(2)) - U_ij*B ).^2*J;
            
            
           DBu=DNu*Nv; DBu=DBu(:);
           DBv=Nu*DNv; DBv=DBv(:);

           DB=[DBu,DBv]/DF;


           H1_err = H1_err + sum((u_Grad(F(1),F(2)) - U_ij*DB ).^2)*J;
           
           
           
D2_B_uu= D2Nu * Nv ; D2_B_uu = D2_B_uu(:);
D2_B_vv= Nu * D2Nv ; D2_B_vv = D2_B_vv(:);
D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:); 



Laplacian_B=D2_B_uu*sum(inv_DF(1,:).^2) + D2_B_uv*2*(inv_DF(1,1)*inv_DF(2,1)+inv_DF(1,2)*inv_DF(2,2))+ ...
            D2_B_vv*sum(inv_DF(2,:).^2) + ...
            DBu*(D_inv_DF_x(1,1) + D_inv_DF_y(1,2)) + DBv*(D_inv_DF_x(2,1) + D_inv_DF_y(2,2));
        
D2_u_h_x_y=D2_B_uu*inv_DF(1,1)*inv_DF(1,2) + D2_B_uv*(inv_DF(1,1)*inv_DF(2,2)+inv_DF(1,2)*inv_DF(2,1))+ ...
            D2_B_vv*inv_DF(2,1)*inv_DF(2,2) + ...
            DBu*D_inv_DF_x(1,2) + DBv*D_inv_DF_x(2,2);


H2_err = H2_err + (Laplacian_u(F(1),F(2)) + D2_u_xy(F(1),F(2)) - U_ij*(Laplacian_B + D2_u_h_x_y)).^2*J;

           

        end
    end
    
    
    


            
end


    L2_err  = sqrt(L2_err );
    H1_err  = sqrt(H1_err );
    H2_err  = sqrt(H2_err );
    
    err = zeros(1,3);
    
    err(1) = L2_err;
    err(2) = H1_err;
    err (3) = H2_err;
    
end


