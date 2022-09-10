function [err,n_dofs]=IGA_2D_Biharmonic(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,test_case)
 
%=====================
% Input:
% ConPts are the control points;
% weights are the weights;
% knotU contains the knot vector in the u direction;
% pu is the degree of NURBS surface in the u direction;
% knotV contains the knot vector in the v direction;
% pv is the degree of NURBS surface in the v direction;
% Refinement mens the times of the h-refinement of both u and v direction;
% t denotes the order to be elevated with analogous to p-refin, i.e., the
% ultimate degree of NURBS  basis functions is (pu+t) in u direction, and
% (pv+t) in the v direction;

% The test_case can be a 'square' and the 'quarter annulus' .

%=====================  
addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')

plot_solution = false;



nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.knotU = knotU; 
nurbs_original.knotV  = knotV; 
nurbs_original.pu = pu;
nurbs_original.pv = pv;


switch test_case
 
    case 'rectangle'
        

    u_Exact=@(x,y,z) exp(x).*exp(y);% Exact solution of Biharmonic equation;
    
    f=@(x,y)  4*exp(x).*exp(y);          % The right hand side of Biharmonic equation;
 
 
    g1=@(x,y) exp(x).*exp(y); % $u|_{\partial \Omega}  = g1$;

    g2=@(x,y) -exp(x).*exp(y); % 在y=0的边界上，梯度u点乘n等于g2;


    u_grad=@(x,y)[exp(x).*exp(y), exp(x).*exp(y)];
                               
     Laplacian_u=@(x,y) 2*exp(x).*exp(y);

     D2_u_xy = @(x,y) exp(x).*exp(y);
    
   

    case 'quarter'
   disp('The domain is a quarter annulus')

    
  u_Exact=@(x,y) exp(x).*exp(y);% Exact solution of Biharmonic equation;
    
  f=@(x,y)  4*exp(x).*exp(y);          % The right hand side of Biharmonic equation;
 
 
  g1=@(x,y) exp(x+y); % 在 r=sqrt(2)/2 处，u =g1;

  g2=@(x,y) -(x+y).*exp(x+y)./sqrt(x.^2+y.^2); % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;



   u_grad=@(x,y)[exp(x).*exp(y), exp(x).*exp(y)];
                               
   Laplacian_u=@(x,y) 2*exp(x).*exp(y);

   D2_u_xy = @(x,y) exp(x).*exp(y);



    case 'disk'
  
    disp('The domain is a unit disk')
     
    u_Exact=@(x,y) sin(x).*sin(y);% Exact solution of Biharmonic equation;
    
    f=@(x,y)  4.*sin(x).*sin(y);          % The right hand side of Biharmonic equation;
 
 
    g1=@(x,y)  sin(x).*sin(y); % 在 r=sqrt(2)/2 处，u =g1;

    g2=@(x,y) (x.*cos(x).*sin(y)+y.*sin(x).*cos(y))./sqrt(x.^2+y.^2); % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;



     u_grad=@(x,y)[cos(x).*sin(y), sin(x).*cos(y)];
                               
     Laplacian_u=@(x,y) -2* sin(x).*sin(y);

     D2_u_xy = @(x,y) cos(x).*cos(y);
   
   

    


    case 'plate_hole'
  
    disp('The domain is a plate hole')
    

     
    u_Exact=@(x,y) sin(x+1).*sin(y+1);% Exact solution of Biharmonic equation;
    
    f=@(x,y)  4.*sin(x+1).*sin(y+1);          % The right hand side of Biharmonic equation;
 
 
    g1=@(x,y)  sin(x+1).*sin(y+1); % 在 r=sqrt(2)/2 处，u =g1;

    g2=@(x,y) (x.*cos(x+1).*sin(y+1)+y.*sin(x+1).*cos(y+1))./sqrt(x.^2+y.^2); % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;



     u_grad=@(x,y)[cos(x+1).*sin(y+1), sin(x+1).*cos(y+1)];
                               
     Laplacian_u=@(x,y) -2* sin(x+1).*sin(y+1);

     D2_u_xy = @(x,y) cos(x+1).*cos(y+1);
     
  
end





if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end


% The following function is used to obtain the information like the global index and the
% coordinates of each element ...
nurbs_refine=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);


disp('=========================================')
disp('The time needed to obtain g1_h and g2_h is:')
tic

[g1_h, g1_bnd_dofs]= L2_project2dirichlet_bnd_g1(ConPts,weights,knotU,pu,knotV,pv,g1,nurbs_refine);

[g2_h, g2_bnd_dofs] =biharmonic_L2_project2dirichlet_bnd_g2(ConPts,weights,knotU,pu,knotV,pv,u_grad,g1_h,nurbs_refine);

toc

n_dofs=nurbs_refine.dof;

u_d_h=zeros(n_dofs,1);

u_d_h(g1_bnd_dofs) = g1_h(g1_bnd_dofs);
u_d_h(g2_bnd_dofs) = g2_h(g2_bnd_dofs);


m=nurbs_refine.m;
n=nurbs_refine.n;

NoEs=nurbs_refine.NoEs;







disp('The degree of  B-spline basis functions is ')
disp(pu)

disp('The number of elements is')
disp(NoEs)

disp('The number of degress of freedom is')
disp(n_dofs)

disp('=========================================')






 [A,rhs] = Solve_2D_Biharmonic_A_F(nurbs_original,nurbs_refine,f); 
 
 


rhs = rhs - A*u_d_h;


[A,rhs]=Iga_2d_bc_biharmonic(A,rhs,m,n);





Uh=A\rhs;

Uh =Uh +u_d_h;




err =Compute_H2_Error(nurbs_original,nurbs_refine,Uh,u_Exact, u_grad,Laplacian_u,D2_u_xy);
 
if plot_solution
plot_uh_2D(nurbs_original, nurbs_refine,Uh,u_Exact)
end

    

end




