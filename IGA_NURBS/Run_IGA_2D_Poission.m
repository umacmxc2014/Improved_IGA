 test_case = 'rectangle';

% test_case = 'quarter';

if strcmp(test_case,'rectangle')

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

t=0;

for i=1:size(Refinement,2)
tStart = clock;
    [err(i,:),dof(i)]=Iga2d_poission(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,test_case);
tEnd = clock;

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

 t=0;
 

% Refinement=[1,2,3,4,5,6,7];
  Refinement=7*ones(1,5);
n_refine=length(Refinement);




for i=1: n_refine
    tic
    ConPts=zeros(3,3,2);
    r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
 knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
pu=2;pv=2;
    [err(i,:),dof(i)]=Iga2d_poission(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,test_case);
toc

end
format short
disp('The CPU time is')
disp(Time)
disp('The average CPU time is')
disp(sum(Time)/length(Time))
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



