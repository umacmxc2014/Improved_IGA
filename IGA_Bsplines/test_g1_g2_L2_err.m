 test_case = 'rectangle';

% test_case = 'quarter';

format short e

if strcmp(test_case,'rectangle')

t=0;  % Modify the degree of the Basis functions.

pu=1;pv=1;



% The ultimate degree is  (pu + t) in the u direction.

ConPts=zeros(2,2,2);
ConPts(:,:,1)=[0 0;1 1];ConPts(:,:,2)=[0 1;0 1];
weights=[1 1;1 1];
knotU=[0 0   1 1];knotV=[0 0  1 1];

Refinement=[2,3,4,5,6];

n_refine=length(Refinement);
err=zeros(n_refine,3);
dof=zeros(n_refine,1);

g1=@(x,y) sin(pi*x)*sin(pi*x);
g2=@(x,y) 4*sin(pi*x)*sin(pi*x);

err=zeros(n_refine,2);

for  i=1:n_refine
[u_h_g1,u_h_g2,err(i,:)] = biharmonic_L2_project2dirichlet_bnd(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,g1,g2);
end

 err_g1=err(:,1); err_g2=err(:,2);
disp('-----------------------')
 disp('The L2 errors for g1 and g2 are')
 disp(err)
 disp('-----------------------')
 
 size(err)
 
 disp('The convergence rate for error in g1 is')
 log(err_g1(1:end-1)./err_g1(2:end))/log(2)
 
 disp('disp(err_g2)=')
 disp(err_g2)

 disp('The convergence rate for error in g2 is')
 
 log(err_g2(1:end-1)./err_g2(2:end))/log(2)

end



if strcmp(test_case,'quarter')


t=2;

ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
 knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
pu=2;pv=2;
Refinement=[2,3,4,5,6];

n_refine=length(Refinement);

dof=zeros(n_refine,1);

 g1=@(x,y) x.^2.*y.^2; % 在 r=sqrt(2)/2 处，u =g1;

 g2=@(x,y) -4*sqrt(2)*x.^2.*y^2; % 在 r=sqrt(2)/2  的边界上，梯度u点乘n等于g2;

err=zeros(n_refine,4);

for  i=1:n_refine
[u_h_g1,u_h_g2,err(i,:)] = biharmonic_L2_project2dirichlet_bnd(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,g1,g2);
end

 err_g1=err(:,1); err_g2=err(:,2);

 disp('The L2 errors for g1 and g2 are')
 disp(err)
 
 disp('The convergence rate for error in g1 is')
 log(err_g1(1:end-1)./err_g1(2:end))/log(2)

 disp('The convergence rate for error in g2 is')
 log(err_g2(1:end-1)./err_g2(2:end))/log(2)


end