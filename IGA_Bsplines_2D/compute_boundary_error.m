function err = compute_boundary_error(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,test_case)




addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')




switch test_case
 
    case 'rectangle'
        

 
 if 1
    g1=@(x,y) exp(x).*exp(y); % $u|_{\partial \Omega}  = g1$;

  
    u_grad=@(x,y)[exp(x).*exp(y), exp(x).*exp(y)];
 end
 
 if 0
     
    g1=@(x,y) sin(pi*x).*(y-2); % $u|_{\partial \Omega}  = g1$;

  
    u_grad=@(x,y)[pi*cos(pi*x)*(y-2), sin(pi*x)];
 end
    
    
   

    case 'quarter'
   disp('The domain is a quarter annulus')

    
 
    g1=@(x,y) exp(x).*exp(y); % $u|_{\partial \Omega}  = g1$;

  
    u_grad=@(x,y)[exp(x).*exp(y), exp(x).*exp(y)];



    case 'disk'
  
    disp('The domain is a unit disk')
     
 
    g1=@(x,y) sin(x).*sin(y); % $u|_{\partial \Omega}  = g1$;

  
    u_grad=@(x,y)[cos(x).*sin(y), sin(x).*cos(y)]; 
   
   

    


    case 'plate_hole'
  
    disp('The domain is a plate hole')
    

     
 
    g1=@(x,y) sin(x+1).*sin(y+1); % $u|_{\partial \Omega}  = g1$;

  
    u_grad=@(x,y) [cos(x+1).*sin(y+1), sin(x+1).*cos(y+1)];
     
  
end





if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end


% The following function is used to obtain the information like the global index and the
% coordinates of each element ...
nurbs_refine=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);


[g1_h]= L2_project2dirichlet_bnd_g1(ConPts,weights,knotU,pu,knotV,pv,g1,nurbs_refine);


[g2_h]= biharmonic_L2_project2dirichlet_bnd_g2(ConPts,weights,knotU,pu,knotV,pv,u_grad,g1_h,nurbs_refine);

err_g1= bnd_g1_error(ConPts,weights,knotU,pu,knotV,pv,g1,g1_h,nurbs_refine);
err_g2= bnd_g2_error(ConPts,weights,knotU,pu,knotV,pv,u_grad,g1_h,g2_h,nurbs_refine);

err = [err_g1,err_g2];


end