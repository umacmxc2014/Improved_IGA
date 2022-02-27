function [W,DW, D2W,F,DF,D2F]=NurbsSurfaceDers(P,U,V,w,pu,u,pv,v)
Uders=bspbasisDers(U,pu,u,2);
Nu=Uders(1,:);
DNu=Uders(2,:);
D2Nu=Uders(3,:);

Vders=bspbasisDers(V,pv,v,2);
Nv=Vders(1,:)';
DNv=Vders(2,:)';
D2Nv=Vders(3,:)';

i=findspan(U,pu,u);
j=findspan(V,pv,v);
u_index=i-pu:i;
v_index=j-pv:j;

w_ij=w(u_index,v_index);

P_u=P(u_index,v_index,1);
P_v=P(u_index,v_index,2);

DIM=2;


F=zeros(DIM,1);



W = Nu * w_ij * Nv;
DW=zeros(1,DIM);
DW(1) = DNu*w_ij*Nv;
DW(2) = Nu*w_ij*DNv;

F(1) = Nu * (w_ij.*P_u)*Nv/W;
F(2) = Nu * (w_ij.*P_v)*Nv/W;

D2W = zeros(DIM,DIM);

D2W(1,1)=D2Nu*w_ij*Nv; D2W(1,2) = DNu * w_ij *DNv;
D2W(2,1) = D2W(1,2);   D2W(2,2) = Nu  * w_ij *D2Nv;

W2 = W*W;
W4 = W2*W2;


R_DNu = (DNu*W-Nu*DW(1))/W2;
R_DNv = (DNv*W-Nv*DW(2))/W2;

DF=zeros(DIM,DIM);

DF(1,1) = R_DNu *(P_u.*w_ij)*Nv;
DF(1,2) = Nu *(P_u.*w_ij)*R_DNv;

DF(2,1) = R_DNu *(P_v.*w_ij)*Nv;
DF(2,2) = Nu *(P_v.*w_ij)*R_DNv;

D2F=zeros(DIM,DIM,DIM);
% D2F(DIM,DIM,1)存储的是 DF关于u的偏导数.
% D2F(DIM,DIM,2)存储的是 DF关于v的偏导数。


R_D2Nu =  ( (D2Nu*W-Nu*D2W(1,1))*W2-(DNu*W-Nu*DW(1))*2*W*DW(1) )/W4; % 这里还没有把权系数 w_{ij} 放进来。
R_D2Nv =  ( (D2Nv*W-Nv*D2W(2,2))*W2-(DNv*W-Nv*DW(2))*2*W*DW(2) )/W4;
R_D2Nuv=  ( (DNu*DW(2) - Nu*D2W(1,2) )*W2- (DNu*W - Nu*DW(1))*2*W*DW(2) )/W4;

D2F(1,1,1) = R_D2Nu*(w_ij.*P_u)*Nv;
D2F(1,2,1) = R_DNu *(w_ij.*P_u)*DNv +  R_D2Nuv * (w_ij.*P_u)*Nv;
D2F(2,1,1) = R_D2Nu*(w_ij.*P_v)*Nv;
D2F(2,2,1) = R_DNu *(w_ij.*P_v)*DNv +  R_D2Nuv * (w_ij.*P_v)*Nv;

D2F(1,1,2) = D2F(1,2,1);
D2F(1,2,2) = Nu *(w_ij.*P_u)*R_D2Nv;
D2F(2,1,2) = D2F(2,2,1);
D2F(2,2,2) = Nu *(w_ij.*P_v)*R_D2Nv;




end
