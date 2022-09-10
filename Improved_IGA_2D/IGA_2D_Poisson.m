function [err, N_DOFs]=IGA_2D_Poisson(ConPts,weights,knotU,pu,knotV,pv,Refinement,t, test_case)

% 现在利用 NURBS 基函数 中的B-spline 基函数作为有限元空间。

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

addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')

nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.knotU = knotU; 
nurbs_original.knotV  = knotV; 
nurbs_original.pu = pu;
nurbs_original.pv = pv;

plot_solution = true;

homogeneous = true;

switch test_case

case 'rectangle'
    
    if homogeneous
u_Exact=@(x,y) sin(pi*x).*sin(pi*y);% Exact solution of Poisson equation;
f=@(x,y)   2*pi*pi*sin(pi*x).*sin(pi*y); % The right hand side of Poisson equation;
u_d=@(x,y) sin(pi*x).*sin(pi*y);
u_Grad=@(x,y)[ pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y) ];
    else
u_Exact=@(x,y)exp(x).*sin(y+1);% Exact solution of Poisson equation;
f=@(x,y)   0; % The right hand side of Poisson equation;
u_d=@(x,y) exp(x).*sin(y+1);
u_Grad=@(x,y)[ exp(x).*sin(y+1), exp(x).*cos(y+1) ];
end


    case 'quarter'
    
u_Exact=@(x,y)exp(x).*exp(y);% Exact solution of Poisson equation;
f=@(x,y)   -2*exp(x).*exp(y); % The right hand side of Poisson equation;
u_d=@(x,y) exp(x).*exp(y);
u_Grad=@(x,y)[ exp(x).*exp(y), exp(x).*exp(y) ];      
    



     case 'disk'
         
u_Exact=@(x,y)exp(x).*exp(y);% Exact solution of Poisson equation;
f=@(x,y)   -2*exp(x).*exp(y); % The right hand side of Poisson equation;
u_d=@(x,y) exp(x).*exp(y);
u_Grad=@(x,y)[ exp(x).*exp(y), exp(x).*exp(y) ];      
    


    case 'plate_hole'
    


u_Exact=@(x,y) sin(x-1).*sin(y-2);% Exact solution of Poisson equation;
f=@(x,y)   2*sin(x-1).*sin(y-2); % The right hand side of Poisson equation;
u_d=@(x,y) sin(x-1).*sin(y-2);
u_Grad=@(x,y)[ cos(x-1).*sin(y-2), sin(x-1).*cos(y-2) ];



end


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
 ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end

nurbs_refine=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);


m=nurbs_refine.m;
n=nurbs_refine.n;

NoEs=nurbs_refine.NoEs;
N_DOFs=nurbs_refine.dof;





[g1_h, ~]= L2_project2dirichlet_bnd_g1(ConPts,weights,knotU,pu,knotV,pv,u_d,nurbs_refine);

% disp('The L2 projection error  is ')
% disp(u_d_h)





disp('================================')
disp('The degree of B-spline basis function is ')
disp(pu)

disp('The number of elements is')
disp(NoEs)

disp('The number of DOFs is')
disp(N_DOFs)

disp('================================')

[A,rhs] = solve_laplace_A_F(nurbs_original,nurbs_refine,f);


% 首先处理下边界 (i.e., v=0) 是 非齐次 Dirichlet 边界。



rhs = rhs - A*g1_h;



[A,rhs]=Iga_2d_bc(A,rhs,m,n);

% spy(A)

U=A\rhs;

% 这时候先把齐次解$u_0$算出来了。
% 现在把非齐次狄利克雷边界的解加上上面的非齐次解里。

U = U+g1_h;





%%
err = Compute_H1_Error(nurbs_original,nurbs_refine,U,u_Exact, u_Grad);

 %%

 
if plot_solution
plot_uh_2D(nurbs_original, nurbs_refine,U,u_Exact)
end

end

