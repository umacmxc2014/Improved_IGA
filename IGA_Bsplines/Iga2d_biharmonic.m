function [err,n_dof]=Iga2d_biharmonic(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,test_case)


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




 
if   strcmp(test_case,'square')





u_Exact=@(x,y) sin(pi*x).*sin(pi*x).*(y-1).^4;
f=@(x,y) -8*( pi^4*cos(2*pi*x).*(y-1).^4 - 6*pi*pi*cos(2*pi*x).*(y-1).^2 ...
                -3*sin(pi*x).*sin(pi*x) ) ;

g1=@(x,y) sin(pi*x).*sin(pi*x); % 在y=0处，u =g1;

g2=@(x,y) 4*sin(pi*x).*sin(pi*x); % 在y=0的边界上，梯度u点乘n等于g2;


u_Exact_grad=@(x,y)[pi*sin(2*pi*x).*(y-1).^4, ...
                                   sin(pi*x).*sin(pi*x)*4*(y-1).^3 ];
                               
Laplacian_u=@(x,y) 2*pi*pi*cos(2*pi*x)*(y-1)^4 + 12*sin(pi*x)*sin(pi*x)*(y-1)^2;

D2_u_xy = @(x,y) 4*pi*sin(2*pi*x).*(y-1).^3;



end


if   strcmp(test_case,'quarter')
  

 
u_Exact=@(x,y) x.^2.*y.^2.*sin(pi*(x.^2+y.^2-2)).*sin(pi*(x.^2+y.^2-2));


  f=@(x,y)4*(256*pi^2*x^2*y^2 + 8*pi^2*(x^4+y^4)-32*pi^4*x^2*y^2*(x^2+y^2)^2-1) *cos(2*pi*(x^2+y^2-2))+...
  64*pi*(x^2+y^2)*(1-12*pi^2*x^2*y^2)*sin(2*pi*(x^2+y^2-2))  + 4;

 g1=@(x,y) x.^2.*y.^2.*sin(pi*(x.^2+y.^2-2)).*sin(pi*(x.^2+y.^2-2)); % 在 r=sqrt(2)/2 处，u =g1;

 g2=@(x,y) -4*sqrt(2)*x.^2.*y^2; % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;


 u_Exact_grad=@(x,y)[y^2*(  x*(1-cos(2*pi*(x^2+y^2-2))) + 2*pi*x^3*sin(2*pi*(x^2+y^2-2)) ), ...
                     x^2*(  y*(1-cos(2*pi*(x^2+y^2-2))) + 2*pi*y^3*sin(2*pi*(x^2+y^2-2)) )];

 Laplacian_u=@(x,y)  y^2*(  1+(8*pi^2*x^4-1)*cos(2*pi*(x^2+y^2-2)) + 10*pi*x^2*sin(2*pi*(x^2+y^2-2))  ) +...
                     x^2*(  1+(8*pi^2*y^4-1)*cos(2*pi*(x^2+y^2-2)) + 10*pi*y^2*sin(2*pi*(x^2+y^2-2))  );

 D2_u_xy = @(x,y) 2*x*y +4*pi*x*y*(x^2+y^2)*sin(2*pi*(x^2+y^2-2)) + 2*x*y*(4*pi^2*x^2*y^2-1)*cos(2*pi*(x^2+y^2-2));


end

%% Project the boundary functions $g_1$ and $g_2$ to the one-dimensional B-splines space.
[u_h_g1,u_h_g2,err] = biharmonic_L2_project2dirichlet_bnd(ConPts,weights,knotU,pu,knotV,pv,Refinement,t,g1,g2);

% disp('error_proj=')
% disp(err)
%%


if t>=1% if t>=1, the degree of NURBS basis fucntions is elevated to (pu+t) in u direction and to (pv+t) in the v direction;
    
[Q,wbar,Ubar,Vbar]=IGADegreeElevSurface(ConPts,weights,knotU,pu,knotV,pv,t);
ConPts=Q;weights=wbar;knotU=Ubar;knotV=Vbar;
 pu=pu+t;pv=pv+t; 
end


% The following function is used to obtain the information like the global index and the
% coordinates of each element ...
nurbsInfo=Iga_2d_grid(knotU,pu,knotV,pv,weights,Refinement);


Element=nurbsInfo.Element;
Coordinate=nurbsInfo.Coordinate;
knotSpanIndex=nurbsInfo.knotSpanIndex;
Ubar=nurbsInfo.Ubar;
Vbar=nurbsInfo.Vbar;
m=nurbsInfo.m;
n=nurbsInfo.n;
Qw=nurbsInfo.Qw;
NoEs=nurbsInfo.NoEs;
n_dof=nurbsInfo.dof;


uNoEs=nurbsInfo.uNoEs;
vNoEs=nurbsInfo.vNoEs;


A=sparse(n_dof,n_dof);rhs=zeros(n_dof,1);
Eledof=(pu+1)*(pv+1);% The number of dofs in an "Element";

np=pu + 1 ;% The number of Gauss quadrature points in the u- or v-parametric direction;



disp('The degree of  B-spline basis functions is ')
disp(pu)

disp('The number of elements is')
disp(NoEs)

disp('The number of degress of freedom is')
disp(n_dof)

disp('=========================================')
tic

for e=1:NoEs % The loop for each element;
    row=Element(e,:); % The global indices of the e-th element.
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);
    uspan=knotSpanIndex(e,1);
   
    vspan=knotSpanIndex(e,2);
    m_row = Eledof;
    n_column = Eledof +1;

	x1={Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv,Eledof};
	x2={np,ue(1),ue(2),ve(1),ve(2),m_row, n_column};
    
    s=Gauss_2d(@(u,v)quad_Ae_Fe(u,v,x1{:},f),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)
  
	Ae=s(1:Eledof,1:Eledof);
    Fe=s(:,end);
	A(row,row)=A(row,row)+Ae;
    rhs(row)=rhs(row)+Fe;
    

    
end

toc

disp('=========================================')

for e=1:(2*uNoEs) % Loop for the  elements along the boundary v=0; 这是双调和的情形
    row=Element(e,:);
    ue=Coordinate(e,1:2);
    ve=Coordinate(e,3:4);
    uspan=knotSpanIndex(e,1);
    vspan=knotSpanIndex(e,2);
    m_row = Eledof;
    n_column =1;
	x1={Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv,Eledof,u_h_g1, u_h_g2};
	x2={np,ue(1),ue(2),ve(1),ve(2), m_row, n_column};
    
    s=Gauss_2d(@(u,v)quad_Fe_bnd(u,v,x1{:} ),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)

    rhs(row)=rhs(row) - s;
end


[A,rhs]=Iga_2d_bc_biharmonic(A,rhs,m,n);


U=A\rhs;

U(1:2*m)=U(1:2*m)+[u_h_g1;u_h_g2];

U=reshape(U,m,n);



  err=[0,0,0];
	    
   for e=1:NoEs
    ue=Coordinate(e,1:2);ve=Coordinate(e,3:4);
    uspan=knotSpanIndex(e,1);vspan=knotSpanIndex(e,2);
 	m_row = 1;
    n_column = 3;
    x1={U,Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv};
 	x2={np,ue(1),ue(2),ve(1),ve(2),m_row, n_column};
   err=err+Gauss_2d(@(u,v)Err_L2(u,v,x1{:},u_Exact,u_Exact_grad,Laplacian_u,D2_u_xy),x2{:});%----s=Gauss_2d(f,np,a,b,c,d)
   end
  	err=sqrt(err);


end


function s=quad_Ae_Fe(u,v,Qw,Ubar,uspan,Vbar,vspan,ConPts,wights,knotU,pu,knotV,pv,Eledof,f)
% ----- Eledof=(pu+1)*(pv+1);
DIM=2;


[W,DW, D2W,F,DF,D2F]=NurbsSurfaceDers(ConPts,knotU,knotV,wights,pu,u,pv,v);
I_2 = eye(2);% 2 by 2 identity matrix.

inv_DF=DF\I_2;

D_DF_x = zeros(DIM,DIM);
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



Uders=bspbasisDers(Ubar,pu,u,2);Vders=bspbasisDers(Vbar,pv,v,2);

Nu=Uders(1,:)';DNu=Uders(2,:)'; D2Nu=Uders(3,:)';
Nv=Vders(1,:) ;DNv=Vders(2,:) ; D2Nv=Vders(3,:); 

%--- i=uspan;j=vspan;

J=abs(det(DF));

B = Nu*Nv;  B=B(:);
s_Fe=f(F(1),F(2))*B*J;


DBu=DNu*Nv; DBu=DBu(:);
DBv=Nu*DNv; DBv=DBv(:);

D2_B_uu= D2Nu * Nv ; D2_B_uu = D2_B_uu(:);
D2_B_vv= Nu * D2Nv ; D2_B_vv = D2_B_vv(:);
D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:); 



Laplacian_B=D2_B_uu*sum(inv_DF(1,:).^2) + D2_B_uv*2*(inv_DF(1,1)*inv_DF(2,1)+inv_DF(1,2)*inv_DF(2,2))+ ...
            D2_B_vv*sum(inv_DF(2,:).^2) + ...
            DBu*(D_inv_DF_x(1,1) + D_inv_DF_y(1,2)) + DBv*(D_inv_DF_x(2,1) + D_inv_DF_y(2,2));



 s_Ae=Laplacian_B*Laplacian_B'*J; % This is the element stiffness matrix for the  biharmonic equation.
  


s=[s_Ae,s_Fe];

end




function s=quad_Fe_bnd(u,v,Qw,Ubar,uspan,Vbar,vspan,ConPts,weights,knotU,pu,knotV,pv,Eledof,u_h_g1, u_h_g2)
% ----- Eledof=(pu+1)*(pv+1);
DIM=2;


[W,DW, D2W,F,DF,D2F]=NurbsSurfaceDers(ConPts,knotU,knotV,weights,pu,u,pv,v);
I_2 = eye(2);% 2 by 2 identity matrix.

inv_DF=DF\I_2;

D_DF_x = zeros(DIM,DIM);
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



Uders=bspbasisDers(Ubar,pu,u,2);Vders=bspbasisDers(Vbar,pv,v,2);

Nu=Uders(1,:)';DNu=Uders(2,:)'; D2Nu=Uders(3,:)';
Nv=Vders(1,:) ;DNv=Vders(2,:) ; D2Nv=Vders(3,:); 
%--- i=uspan;j=vspan;
J=abs(det(DF));

B = Nu*Nv;  B=B(:);



DBu=DNu*Nv; DBu=DBu(:);
DBv=Nu*DNv; DBv=DBv(:);

D2_B_uu= D2Nu * Nv ; D2_B_uu = D2_B_uu(:);
D2_B_vv= Nu * D2Nv ; D2_B_vv = D2_B_vv(:);
D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:); 



Laplacian_B=D2_B_uu*sum(inv_DF(1,:).^2) + D2_B_uv*2*(inv_DF(1,1)*inv_DF(2,1)+inv_DF(1,2)*inv_DF(2,2))+ ...
            D2_B_vv*sum(inv_DF(2,:).^2) + ...
            DBu*(D_inv_DF_x(1,1) + D_inv_DF_y(1,2)) + DBv*(D_inv_DF_x(2,1) + D_inv_DF_y(2,2));


 u_dof = pu+1;

 if vspan==pv+1 % 这说明当前的单元在边界的第一层单元上
 s_Ae=Laplacian_B*( Laplacian_B(1:2*u_dof)' * [u_h_g1(uspan-pu:uspan); u_h_g2(uspan-pu:uspan)])*J; % This is the element stiffness matrix for the  biharmonic equation.
 
 else if vspan==pv+2 % 这说明当前的单元在边界的第二层单元上
 
 s_Ae=Laplacian_B*( [zeros(1,u_dof),Laplacian_B(1:u_dof)'] * [u_h_g1(uspan-pu:uspan); u_h_g2(uspan-pu:uspan)])*J; % This is the element stiffness matrix for the  biharmonic equation.
     else
         disp('The wrong elements have been used!')
  
     end
 end


s=s_Ae;

end



function s=Err_L2(u,v,U,Qw,Ubar,uspan,Vbar,vspan,ConPts,wights,knotU,pu,knotV,pv,u_Exact,u_Exact_grad,Laplacian_u,D2_u_xy)

% ----- Eledof=(pu+1)*(pv+1);
DIM=2;

Eledof=(pu+1)*(pv+1);

[W,DW, D2W,F,DF,D2F]=NurbsSurfaceDers(ConPts,knotU,knotV,wights,pu,u,pv,v);
I_2 = eye(2);% 2 by 2 identity matrix.

inv_DF=DF\I_2;

D_DF_x = zeros(DIM,DIM);
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



Uders=bspbasisDers(Ubar,pu,u,2);Vders=bspbasisDers(Vbar,pv,v,2);

Nu=Uders(1,:)';DNu=Uders(2,:)'; D2Nu=Uders(3,:)';
Nv=Vders(1,:) ;DNv=Vders(2,:) ; D2Nv=Vders(3,:); 
%--- i=uspan;j=vspan;

J=abs(det(DF));
U_ij = U(uspan-pu:uspan,vspan-pv:vspan);
U_ij = U_ij(:);
U_ij = U_ij'; % 列向量

B = Nu*Nv;  B=B(:);

s_L2_err = (u_Exact(F(1),F(2)) - U_ij*B  ).^2*J;

DBu=DNu*Nv; DBu=DBu(:);
DBv=Nu*DNv; DBv=DBv(:);

DB=[DBu,DBv]/DF;


s_H1_err = sum((u_Exact_grad(F(1),F(2)) - U_ij*DB ).^2)*J;


D2_B_uu= D2Nu * Nv ; D2_B_uu = D2_B_uu(:);
D2_B_vv= Nu * D2Nv ; D2_B_vv = D2_B_vv(:);
D2_B_uv= DNu*DNv;    D2_B_uv = D2_B_uv(:); 



Laplacian_B=D2_B_uu*sum(inv_DF(1,:).^2) + D2_B_uv*2*(inv_DF(1,1)*inv_DF(2,1)+inv_DF(1,2)*inv_DF(2,2))+ ...
            D2_B_vv*sum(inv_DF(2,:).^2) + ...
            DBu*(D_inv_DF_x(1,1) + D_inv_DF_y(1,2)) + DBv*(D_inv_DF_x(2,1) + D_inv_DF_y(2,2));
        
D2_u_h_x_y=D2_B_uu*inv_DF(1,1)*inv_DF(1,2) + D2_B_uv*(inv_DF(1,1)*inv_DF(2,2)+inv_DF(1,2)*inv_DF(2,1))+ ...
            D2_B_vv*inv_DF(2,1)*inv_DF(2,2) + ...
            DBu*D_inv_DF_x(1,2) + DBv*D_inv_DF_x(2,2);


s_H2_err = (Laplacian_u(F(1),F(2)) + D2_u_xy(F(1),F(2)) - U_ij*(Laplacian_B + D2_u_h_x_y)).^2*J;


s=[s_L2_err,s_H1_err,s_H2_err];
  



end

% Example:
%=========================================================
% t=0;
% knotU=[0,0,0,1,1,1];
% knotV=knotU;
% a=sqrt(2)/2;
% wights=[1,a,1;a,1,a;1,a,1];
% ConPts=zeros(3,3,2);
% ConPts(:,:,1)=[-a,0,a; -2*a   0   2*a;-a, 0, a];
% ConPts(:,:,2)=[ a,2*a,  a; 0,0, 0;-a, -2*a, -a];
% pu=2;pv=2;
% Refinement=3;
% tic
% [err,dof]=Iga2d_biharmonic(ConPts,wights,knotU,pu,knotV,pv,Refinement,t)
% toc
