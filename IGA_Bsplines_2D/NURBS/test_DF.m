u=0:0.005:1;
v=u;
[u,v] = meshgrid(u,v);

%%

pu=2;pv=2;

t= 1 ;

% Elevate the degree of the B-plines basis functions to (pu+t) and (pv+t)
% For the biharmonic equation, the t should be $t>=0$ for the quarter annulus domain.

ConPts=zeros(3,3,2);
r=sqrt(2)/2; R=sqrt(2); rR=(r+R)/2;

ConPts(:,:,1)=[0 0 0;r rR R;r rR  R]; 
ConPts(:,:,2)=[r rR R;r rR R;0 0 0];
weights=[2 2 2;1 1 1;1 1 1];
knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
%%

%%
pu=2;pv=2;

t= 3;

% Elevate the degree of the B-plines basis functions to (pu+t) and (pv+t)
% For the biharmonic equation, the t should be $t>=0$ for the quarter annulus domain.

ConPts=zeros(3,3,2);
a=sqrt(2)/2;

ConPts(:,:,1)=[-a      0  a; -2*a   0   2*a; -a       0   a]; 
ConPts(:,:,2)=[ a   2*a  a;      0   0      0; -a  -2*a  -a];
weights=[1,a,1;a,1,a;1,a,1];
knotU=[0 0  0 1 1 1];knotV=[0 0 0 1 1 1];
%%
controlPts        = zeros(4,3,2);

controlPts(1,1,:) = [0 1];
controlPts(1,2,:) = [0 2.5];
controlPts(1,3,:) = [0 4];

controlPts(2,1,:) = [0.4142135623730951 1];
controlPts(2,2,:) = [1.5 2.5];
controlPts(2,3,:) = [4 4];

controlPts(3,1,:) =  [1 0.4142135623730951];
controlPts(3,2,:) = [2.5  1.5];
controlPts(3,3,:) = [4 4];

controlPts(4,1,:) = [1 0];
controlPts(4,2,:) = [2.5 0];
controlPts(4,3,:) = [4 0];

ConPts=controlPts;

knotU=[0 0 0  0.5 1 1 1];pu=2;
knotV=[0 0 0  1 1 1];  pv=2;

a= 0.5*(1+1/sqrt(2));
weights=[1 1 1;a a 1;a a 1;1 1 1];

%%



[m,n] =size(u);
det_DF=zeros(m,n);

for i=1:m
    for j=1:n
         [S,DF,W,DWu,DWv]=NurbsSurface(ConPts,weights,knotU,pu,u(i,j),knotV,pv,v(i,j));
         det_DF(i,j) = det(DF);
         if(abs(det_DF(i,j))<=1.0e-5)
             disp('=======')
             disp('u and v are')
             disp(u(i,j))
             disp(v(i,j))
             disp(S)
             disp('=======')
         end
    end
end

surf(u,v,det_DF)
colorbar