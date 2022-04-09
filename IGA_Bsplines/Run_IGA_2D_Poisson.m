 test_case = 'rectangle';

 % test_case = 'quarter';

if strcmp(test_case,'rectangle')
    
t=0;  % The (pu+t) is the ultimate degree of B splines basis functions
%% Case I: rectangle domain
ConPts=zeros(2,2,2);
ConPts(:,:,1)=[0 0;1 1];ConPts(:,:,2)=[0 1;0 1];
weights=[1 1;1 1];
 knotU=[0 0   1 1];knotV=[0 0  1 1];
pu=1;pv=1;
Refinement=[1,2,3,4,5];

n_refine=length(Refinement);
err=zeros(n_refine,2);
dof=zeros(n_refine,1);



for i=1:size(Refinement,2)
tic
    [err(i,:),dof(i)]=Iga2d_poisson(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,test_case);
    
toc
end
format short e
disp('The degree of NURBS basis is ')
disp(pu+t)
 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error'  ])
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(dof)
disp('-----------------------------------------------------------------')
  disp('The convergence order in L2, H1 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2)])
 disp('================================================================')



%%

end

%% Case I: quarter domain

if strcmp(test_case,'quarter')

ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
pu=2;pv=2;


t=0;  % Elevate the degree of B-splines to (pu+t).
    




Refinement=[1,2,3,4,5,6];

% Refinement=[7,8];

n_refine=length(Refinement);

err=zeros(n_refine,2);
dof=zeros(n_refine,1);
Time=zeros(n_refine,1);


for i=1: n_refine
tic
    [err(i,:),dof(i)]=Iga2d_poisson(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,test_case);
toc
end



format short e
disp('The degree of NURBS basis is ')
disp(pu+t)

 disp('================================================================')
 disp(['The L2 norm error            ||       The H1 norm error'  ])
 
 disp('----------------------------------------------------------------')
 disp([err])
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(dof)
disp('-----------------------------------------------------------------')
  disp('The convergence order in L2, H1 are:  ')
 
  disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2)])
 disp('================================================================')

end



