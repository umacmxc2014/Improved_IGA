function [err,n_dofs]=IGA_2D_Biharmonic_Homogeneous(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,test_case)
 

%%

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

%% 
addpath('./IGA_Grid_data/')
addpath('./NURBS/')
addpath('./quadrature/')

plot_solution = true;



nurbs_original.ConPts = ConPts;
nurbs_original.weights = weights;
nurbs_original.knotU = knotU; 
nurbs_original.knotV  = knotV; 
nurbs_original.pu = pu;
nurbs_original.pv = pv;



 

switch test_case
 
    case 'rectangle'

     pi_4 = pi^4;
    
    u_Exact=@(x,y,z) sin(pi*x).*sin(pi*y).* sin(pi*x).*sin(pi*y);% Exact solution of Biharmonic equation;
    
     f=@(x,y)  -8*pi_4*(cos(2*pi*x)*sin(pi*y)*sin(pi*y)   + ...
                                      cos(2*pi*y)*sin(pi*x)*sin(pi*x)) + ...
                         8*pi_4*cos(2*pi*x)*cos(2*pi*y) ;              % The right hand side of Biharmonic equation;
 
 

    g1=@(x,y) sin(pi*x)*sin(pi*y)* sin(pi*x)*sin(pi*y); % 在y=0处，u =g1;

    g2=@(x,y) 0; % 在y=0的边界上，梯度u点乘n等于g2;


    u_grad=@(x,y)[pi*sin(2*pi*x)*sin(pi*y)*sin(pi*y) , ...
                            pi*sin(2*pi*y)*sin(pi*x)*sin(pi*x)];
                               
     Laplacian_u=@(x,y) 2*pi*pi*( cos(2*pi*x)*sin(pi*y)*sin(pi*y) + cos(2*pi*y)*sin(pi*x)*sin(pi*x) );

     D2_u_xy = @(x,y) pi*pi*sin(2*pi*x)*sin(2*pi*y);

    
   

    case 'quarter'
   disp('The domain is a quarter annulus')

  u_Exact=@(x,y) x.^2.*y.^2.*sin(pi*(x.^2+y.^2-2)).*sin(pi*(x.^2+y.^2-2));

  f=@(x,y)4*(256*pi^2*x^2*y^2 + 8*pi^2*(x^4+y^4)-32*pi^4*x^2*y^2*(x^2+y^2)^2-1) *cos(2*pi*(x^2+y^2-2))+...
  64*pi*(x^2+y^2)*(1-12*pi^2*x^2*y^2)*sin(2*pi*(x^2+y^2-2))  + 4;

 g1=@(x,y) x.^2.*y.^2.*sin(pi*(x.^2+y.^2-2)).*sin(pi*(x.^2+y.^2-2)); % 在 r=sqrt(2)/2 处，u =g1;

g2=@(x,y) -4*sqrt(2)*x.^2.*y^2; % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;


 u_grad=@(x,y)[y^2*(  x*(1-cos(2*pi*(x^2+y^2-2))) + 2*pi*x^3*sin(2*pi*(x^2+y^2-2)) ), ...
                     x^2*(  y*(1-cos(2*pi*(x^2+y^2-2))) + 2*pi*y^3*sin(2*pi*(x^2+y^2-2)) )];

 Laplacian_u=@(x,y)  y^2*(  1+(8*pi^2*x^4-1)*cos(2*pi*(x^2+y^2-2)) + 10*pi*x^2*sin(2*pi*(x^2+y^2-2))  ) +...
                     x^2*(  1+(8*pi^2*y^4-1)*cos(2*pi*(x^2+y^2-2)) + 10*pi*y^2*sin(2*pi*(x^2+y^2-2))  );

 D2_u_xy = @(x,y) 2*x*y +4*pi*x*y*(x^2+y^2)*sin(2*pi*(x^2+y^2-2)) + 2*x*y*(4*pi^2*x^2*y^2-1)*cos(2*pi*(x^2+y^2-2));




    case 'disk'
  
    disp('The domain is a disk')
  
     
    u_Exact=@(x,y) 1 - cos(2*(x.^2+y.^2-1));% Exact solution of Biharmonic equation;
  
    % The right hand side of Biharmonic equation;  
    f=@(x,y)  16*(  (3-16*x.^4).*cos(2*(x.^2+y.^2-1)) - 24*x.^2*sin(2*(x.^2+y.^2-1)) )+...
                    16*(  (3-16*y.^4).*cos(2*(x.^2+y.^2-1)) - 24*y.^2*sin(2*(x.^2+y.^2-1)) )+... ;          
                    2*16*((1-16*x.^2.*y.^2).*cos(2*(x.^2+y.^2-1))-4*(x.^2+y.^2).*sin(2*(x.^2+y.^2-1)));
 
    g1=@(x,y) 0; % 在 r=sqrt(2)/2 处，u =g1;

    g2=@(x,y) 0; % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;



    u_grad=@(x,y)4*[sin(2*(x.^2+y.^2-1)).*x, sin(2*(x.^2+y.^2-1)).*y];
                               
    Laplacian_u=@(x,y) 4*(sin(2*(x.^2+y.^2-1))  + 4*x.^2.*cos(2*(x.^2+y.^2-1)) )+...
                                    4*(sin(2*(x.^2+y.^2-1))  + 4*y.^2.*cos(2*(x.^2+y.^2-1)) ) ;

     D2_u_xy = @(x,y) 16*x.*y.*cos(2*(x.^2+y.^2-1));    
     
     
     u_Exact = @(x,y) (x.^2+y.^2-1).^4;
     f=@(x,y) 1152*x^2*(x^2 + y^2 - 1) + 144*(x^2 + y^2 - 1)^2 + 384*x^4 +...
         1152*y^2*(x^2 + y^2 - 1) + 144*(x^2 + y^2 - 1)^2 + 384*y^4 + ...
         +2*(  384*x^2*y^2 + 192*x^2*(x^2 + y^2 - 1) + 192*y^2*(x^2 + y^2 - 1) + 48*(x^2 + y^2 - 1)^2 );
     
     u_grad=@(x,y) [8*x*(x^2 + y^2 - 1)^3,8*y*(x^2 + y^2 - 1)^3];
     Laplacian_u=@(x,y)8*(x^2 + y^2 - 1)^3 + 48*x^2*(x^2 + y^2 - 1)^2 + ...
                                    8*(x^2 + y^2 - 1)^3 + 48*y^2*(x^2 + y^2 - 1)^2;


    D2_u_xy = @(x,y) 48*x*y*(x^2 + y^2 - 1)^2;



    case 'plate_hole'
  
    disp('The domain is a plate hole')
    

     
   % u_Exact=@(x,y) (1-cos(2*pi*x)).*(1-cos(2*pi*y)).*( 1-cos(2*(x.^2+y.^2-1))  );% Exact solution of Biharmonic equation;
    
  %  f=@(x,y) 48*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 16*pi^4*cos(2*pi*x)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*y) - 1) - 256*x^4*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 384*x^2*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 96*pi^2*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*x)*(cos(2*pi*y) - 1) + 128*x*pi^3*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) + 512*x^3*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) - 384*x^2*pi^2*cos(2*x^2 + 2*y^2 - 2)*cos(2*pi*x)*(cos(2*pi*y) - 1) - 384*x*pi*cos(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) + ...
  %      48*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 16*pi^4*cos(2*pi*y)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*x) - 1) - 256*y^4*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 384*y^2*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 96*pi^2*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*y)*(cos(2*pi*x) - 1) + 128*y*pi^3*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1) + 512*y^3*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1) - 384*y^2*pi^2*cos(2*x^2 + 2*y^2 - 2)*cos(2*pi*y)*(cos(2*pi*x) - 1) - 384*y*pi*cos(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1) + ...
  %      2*( 16*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 64*x^2*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 64*y^2*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 16*pi^4*cos(2*pi*x)*cos(2*pi*y)*(cos(2*x^2 + 2*y^2 - 2) - 1) - 16*pi^2*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*x)*(cos(2*pi*y) - 1) - 16*pi^2*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*y)*(cos(2*pi*x) - 1) - 256*x^2*y^2*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 64*x^2*pi^2*cos(2*x^2 + 2*y^2 - 2)*cos(2*pi*y)*(cos(2*pi*x) - 1) - 64*y^2*pi^2*cos(2*x^2 + 2*y^2 - 2)*cos(2*pi*x)*(cos(2*pi*y) - 1) - 64*x*pi*cos(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) - 64*y*pi*cos(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1) + 64*x*pi^3*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*y)*sin(2*pi*x) + 64*y*pi^3*sin(2*x^2 + 2*y^2 - 2)*cos(2*pi*x)*sin(2*pi*y) + 256*x*y*pi^2*cos(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*sin(2*pi*y) + 256*x*y^2*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) + 256*x^2*y*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1)); 


  
  u_Exact = @(x,y) x.^2.*(4-x).^2.*y.^2.*(4-y).^2.*(x.^2+y.^2-1).^2;
  
  f=@(x,y)(48*(2*x - 4)^2*(x^2 + y^2 - 1) + 24*(x^2 + y^2 - 1)^2 - 192*x^2*(- x^2 + 4*x) + 24*(- x^2 + 4*x)^2 + 96*x^2*(2*x - 4)^2 - 96*(- x^2 + 4*x)*(x^2 + y^2 - 1) + 64*x*(2*x - 4)*(x^2 + y^2 - 1) + 32*x*(8*x - 16)*(x^2 + y^2 - 1) - 192*x*(2*x - 4)*(- x^2 + 4*x))*y^2*(4-y)^2 +...
                x^2*(4-x)^2*(48*(2*y - 4)^2*(x^2 + y^2 - 1) + 24*(x^2 + y^2 - 1)^2 - 192*y^2*(- y^2 + 4*y) + 24*(- y^2 + 4*y)^2 + 96*y^2*(2*y - 4)^2 - 96*(- y^2 + 4*y)*(x^2 + y^2 - 1) + 64*y*(2*y - 4)*(x^2 + y^2 - 1) + 32*y*(8*y - 16)*(x^2 + y^2 - 1) - 192*y*(2*y - 4)*(- y^2 + 4*y))+...
                2*(8*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2 - 32*x^2*(- x^2 + 4*x)^2*(- y^2 + 4*y) - 32*y^2*(- x^2 + 4*x)*(- y^2 + 4*y)^2 + 16*x^2*(2*y - 4)^2*(- x^2 + 4*x)^2 + 16*y^2*(2*x - 4)^2*(- y^2 + 4*y)^2 + 16*(- x^2 + 4*x)*(- y^2 + 4*y)*(x^2 + y^2 - 1)^2 - 16*(- x^2 + 4*x)*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) - 16*(- x^2 + 4*x)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1) - 8*(2*y - 4)^2*(- x^2 + 4*x)*(x^2 + y^2 - 1)^2 + 8*(2*y - 4)^2*(- x^2 + 4*x)^2*(x^2 + y^2 - 1) - 8*(2*x - 4)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1)^2 + 8*(2*x - 4)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) + 4*(2*x - 4)^2*(2*y - 4)^2*(x^2 + y^2 - 1)^2 - 32*x*(2*x - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)^2 - 32*y*(2*y - 4)*(- x^2 + 4*x)^2*(- y^2 + 4*y) - 32*x*(2*x - 4)*(2*y - 4)^2*(- x^2 + 4*x)*(x^2 + y^2 - 1) - 32*y*(2*x - 4)^2*(2*y - 4)*(- y^2 + 4*y)*(x^2 + y^2 - 1) + 64*x*(2*x - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)*(x^2 + y^2 - 1) + 64*y*(2*y - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)*(x^2 + y^2 - 1) + 128*x*y*(2*x - 4)*(2*y - 4)*(- x^2 + 4*x)*(- y^2 + 4*y));

        




                 % The right hand side of Biharmonic equation;
 
 
    g1=@(x,y)  0; % 在 r=sqrt(2)/2 处，u =g1;

    g2=@(x,y) 0; % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;



%     u_grad=@(x,y)[2*pi*sin(2*pi*x)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*y) - 1) + 4*x*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1),...
%        2*pi*sin(2*pi*y)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*x) - 1) + 4*y*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1)];

                             

                               
%     Laplacian_u=@(x,y) 4*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) + 4*pi^2*cos(2*pi*x)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*y) - 1) + 16*x^2*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 16*x*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1) + ...
%         4*sin(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) + 4*pi^2*cos(2*pi*y)*(cos(2*x^2 + 2*y^2 - 2) - 1)*(cos(2*pi*x) - 1) + 16*y^2*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 16*y*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1);


%     D2_u_xy = @(x,y)  16*x*y*cos(2*x^2 + 2*y^2 - 2)*(cos(2*pi*x) - 1)*(cos(2*pi*y) - 1) - 4*pi^2*sin(2*pi*x)*sin(2*pi*y)*(cos(2*x^2 + 2*y^2 - 2) - 1) - 8*x*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*y)*(cos(2*pi*x) - 1) - 8*y*pi*sin(2*x^2 + 2*y^2 - 2)*sin(2*pi*x)*(cos(2*pi*y) - 1);

u_grad = @(x,y)[4*x*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) - 2*(2*x - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)^2*(x^2 + y^2 - 1)^2,...
    4*y*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) - 2*(2*y - 4)*(- x^2 + 4*x)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1)^2];
     
Laplacian_u=@(x,y) 4*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) - 4*(- x^2 + 4*x)*(- y^2 + 4*y)^2*(x^2 + y^2 - 1)^2 + 2*(2*x - 4)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1)^2 + 8*x^2*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2 - 16*x*(2*x - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) +... 
    4*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2*(x^2 + y^2 - 1) - 4*(- x^2 + 4*x)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1)^2 + 2*(2*y - 4)^2*(- x^2 + 4*x)^2*(x^2 + y^2 - 1)^2 + 8*y^2*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2 - 16*y*(2*y - 4)*(- x^2 + 4*x)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1);


 D2_u_xy = @(x,y) 8*x*y*(- x^2 + 4*x)^2*(- y^2 + 4*y)^2 + 4*(2*x - 4)*(2*y - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)*(x^2 + y^2 - 1)^2 - 8*x*(2*y - 4)*(- x^2 + 4*x)^2*(- y^2 + 4*y)*(x^2 + y^2 - 1) - 8*y*(2*x - 4)*(- x^2 + 4*x)*(- y^2 + 4*y)^2*(x^2 + y^2 - 1);

end


 


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end





% The following function is used to obtain the information like the global index and the
% coordinates of each element ...
nurbs_refine=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);




 
NoEs = nurbs_refine.NoEs;


n_dofs = nurbs_refine.n_dofs;

m = nurbs_refine.m;

n = nurbs_refine.n;



disp('=========================================')

disp(['The degree of  B-spline basis functions is:  ',num2str(pu)])
disp(['The number of elements is:  ',num2str(NoEs)])

disp(['The number of degress of freedom is: ',num2str(n_dofs)])

disp('=========================================')



tic

 [A,rhs] = Solve_2D_Biharmonic_A_F(nurbs_original,nurbs_refine,f); 
 
 toc


[A,rhs]=Iga_2d_bc_biharmonic(A,rhs,m,n);





Uh=A\rhs;




err =Compute_H2_Error(nurbs_original,nurbs_refine,Uh,u_Exact, u_grad,Laplacian_u,D2_u_xy);
 
if plot_solution
plot_uh_2D(nurbs_original, nurbs_refine,Uh,u_Exact)
end

    

end




