
 % test_case = 'square';

 test_case = 'quarter';


% The test_case can be 'square' and the 'quarter' ;


if strcmp(test_case, 'square')

%% Case I: the domain is an unit square.
% The physical domain for example is $[0,1]^2$.


pu=1;pv=1;

t=1;  % Modify the degree of the Basis functions.

% The ultimate degree is  (pu + t) in the u direction.

ConPts=zeros(2,2,2);
ConPts(:,:,1)=[0 0;1 1];ConPts(:,:,2)=[0 1;0 1];
weights=[1 1;1 1];
knotU=[0 0   1 1];knotV=[0 0  1 1];

Refinement=[1,2,3,4,5];

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
 disp('----------------------------------------------------------------')
 disp('The convergence order in L2, H1 and H2 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2),log(err(1:end-1,3)./err(2:end,3))/log(2)])
 disp('================================================================')

end


if strcmp(test_case, 'quarter' )

%% Case II.
% The physical domain  is a quarter annulus.


t=3 ;

ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
 knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
pu=2;pv=2;

Refinement=[2,3,4,5,6];

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
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(dof)
disp('-----------------------------------------------------------------')
  disp('The convergence order in L2, H1 and H2 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2),log(err(1:end-1,3)./err(2:end,3))/log(2)])
 disp('================================================================')

 end
