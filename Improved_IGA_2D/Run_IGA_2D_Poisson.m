%%  Currently, there are 4 test cases.

%  test_case = 'rectangle';

%  test_case = 'quarter';

% test_case = 'disk';
 
  test_case = 'plate_hole';

%%
t= 2; % Elevate the degree of B-spline basis function to (pu+t) in u-direction, and to (pv+t) in the v-direction
%%

switch test_case

case 'rectangle'  
%% Case I: rectangular domain 
pu=1;pv=1; % The ultimate degree of B splines basis functions used is (pu+t).
ConPts=zeros(2,2,2);
ConPts(:,:,1)=[0 0;1 1];ConPts(:,:,2)=[0 1;0 1];
weights=[1 1;1 1];
knotU=[0 0  1 1];
knotV=[0 0  1 1];

case 'quarter'
%% Case II: quarter domain    
pu=2;pv=2; % The ultimate degree of B splines basis functions used is (pu+t).
ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;
ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; 
ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
knotU=[0 0  0 1 1 1];
knotV=[0 0 0 1 1 1];

%%



case  'disk'
%% Case III: the domain is a unit disk

pu=2;pv=2;  % The ultimate degree of B splines basis functions used is (pu+t).
ConPts=zeros(3,3,2);
a=sqrt(2)/2;
ConPts(:,:,1)=[-a      0  a; -2*a   0   2*a; -a       0   a]; 
ConPts(:,:,2)=[ a   2*a  a;      0   0      0; -a  -2*a  -a];
weights=[1,a,1;a,1,a;1,a,1];
knotU=[0 0  0 1 1 1];
knotV=[0 0 0 1 1 1];

case  'plate_hole'

%% Case IV.
% The physical domain for example is a plate with a hole.

pu=2; pv=2; % The ultimate degree of B splines basis functions used is (pu+t).
ConPts        = zeros(4,3,2);

ConPts(1,1,:) = [0 1];
ConPts(1,2,:) = [0 2.5];
ConPts(1,3,:) = [0 4];

ConPts(2,1,:) = [0.4142135623730951 1];
ConPts(2,2,:) = [1.5 2.5];
ConPts(2,3,:) = [4 4];

ConPts(3,1,:) =  [1 0.4142135623730951];
ConPts(3,2,:) = [2.5  1.5];
ConPts(3,3,:) = [4 4];

ConPts(4,1,:) = [1 0];
ConPts(4,2,:) = [2.5 0];
ConPts(4,3,:) = [4 0];


knotU=[0 0 0  0.5 1 1 1];
knotV=[0 0 0  1 1 1];   

a= 0.5*(1+1/sqrt(2));
weights=[1 1 1;a 0.8 1;a 0.8 1;1 1 1];
end
 
%%

Refinement = [8]; % h-refine the initial mesh 1,2,3,4,5,6,7 times.

n_refine=length(Refinement);

err=zeros(n_refine,2);
dof=zeros(n_refine,1);
Time=zeros(n_refine,1);

for i=1: n_refine
tic
    [err(i,:), dof(i)]=IGA_2D_Poisson(ConPts,weights,knotU,pu,knotV,pv,Refinement(i),t,test_case);
toc
end
  


 disp('================================================================')
 disp('The L2 norm error            ||       The H1 norm error       || The Laplacian u - Laplacian u_h')
 format short e 
 disp('----------------------------------------------------------------')
 disp(err)
 disp('----------------------------------------------------------------')
disp('The #DOFs:')
disp(dof)

format short
disp('-----------------------------------------------------------------')
disp('The convergence order in L2, H1 and H2 are:  ')
 L2_order = log(err(1:end-1,1)./err(2:end,1))/log(2);
 H1_order = log(err(1:end-1,2)./err(2:end,2))/log(2);
 disp([log(err(1:end-1,1)./err(2:end,1))/log(2),log(err(1:end-1,2)./err(2:end,2))/log(2)])
 disp('================================================================')
 
 
 
sprintf('%.3e %s %s%.3e %s %s %s', err(1,1), '[-]', '&',err(1,2), '[-]', '&')
 for i=2:n_refine
     sprintf('%.3e %s%.2f%s %s%.3e %s%.2f%s %s', err(i,1), '[', L2_order(i-1), ']', '&',err(i,2),'[', H1_order(i-1), ']', '&')  
 end
     

%%


