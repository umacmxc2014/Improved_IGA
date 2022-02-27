
% test_case = 'square';

test_case = 'quarter';

% The test_case can be 'square' and the 'quarter' ;

format short e

if strcmp(test_case, 'square')

%% Case I.
% The physical domain a square domain  $\Omega = [0,1]^2$.


pu=1;pv=1;

t=1; % $t \ge 1$ for the polygonal domain. 
% Elevate the degree of the B-plines basis functions.
% For the biharmonic equation, the t should be t>=1 for the polygonal domain.

% The ultimate degree is  (pu + t) in the u direction.

ConPts=zeros(2,2,2);
ConPts(:,:,1)=[0 0;1 1];ConPts(:,:,2)=[0 1;0 1];
weights=[1 1;1 1];
knotU=[0 0   1 1];knotV=[0 0  1 1];

Refinement=[2,3,4,5,6]; % h-refine the initial mesh 2,3,4,5,6 times.

n_refine=length(Refinement);
err=zeros(n_refine,3);
dof=zeros(n_refine,1);


 for  i=1:n_refine
 [err(i,:),dof(i)]=Iga2d_biharmonic(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t, test_case);
 end
 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error       || The Laplacian u - Laplacian u_h'])
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
 disp('The #DOFs:')
 disp(dof)


 format short
 disp('----------------------------------------------------------------')
 disp('The convergence order in L2, H1 and H2 are:  ')

 disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2),log(err(1:end-1,3)./err(2:end,3))/log(2)])
 disp('================================================================')

end


if strcmp(test_case, 'quarter')

%% Case II.
% The physical domain for example is a quarter annulus.

pu=2;pv=2;

t=0;
% Elevate the degree of the B-plines basis functions to (pu+t) and (pv+t)
% For the biharmonic equation, the t should be $t>=0$ for the quarter annulus domain.

ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];

Refinement = [2,3,4,5,6]; % h-refine the initial mesh 2,3,4,5,6 times.

n_refine=length(Refinement);

err=zeros(n_refine,3);
dof=zeros(n_refine,1);
Time=zeros(n_refine,1);

  for  i=1:n_refine
     tic
 [err(i,:),dof(i)]=Iga2d_biharmonic(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t, test_case);
     toc
  end
  


 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error       || The Laplacian u - Laplacian u_h'])
 format short e 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(dof)

format short
disp('-----------------------------------------------------------------')
disp('The convergence order in L2, H1 and H2 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2),log(err(1:end-1,3)./err(2:end,3))/log(2)])
 disp('================================================================')

 end
