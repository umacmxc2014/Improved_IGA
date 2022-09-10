function plot_uh_2D(nurbs_original,nurbs_refine,Uh,u_Exact)
n_pnts=200;
m_pnts=n_pnts;

ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;

knotU=nurbs_refine.Ubar;
knotV=nurbs_refine.Vbar;




pu            =  nurbs_refine.pu;
pv            =  nurbs_refine.pv;

m = length(knotU) - pu -1;
n = length(knotV) - pv -1;

UBreaks=nurbs_refine.UBreaks;  
VBreaks=nurbs_refine.VBreaks;    

u_pnts = linspace(UBreaks(1),UBreaks(end),m_pnts+1);
v_pnts = linspace(VBreaks(1),VBreaks(end),n_pnts+1);

X = zeros(m_pnts+1,n_pnts+1);
Y = X;

 for j=1:(n_pnts+1)
    for i=1:(m_pnts+1)
     S= PointOnNurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u_pnts(i),knotV_o,pv_o,v_pnts(j));
     X(i,j) = S(1);
     Y(i,j) = S(2);
    end
 end

X = X';
Y = Y';


Z = zeros(m_pnts+1,n_pnts+1);

for j=1:(n_pnts+1)
     v_span = findspan(knotV,pv,v_pnts(j));
     Nv = bspbasisDers(knotV,pv,v_pnts(j),0);
     Nv = Nv(1,:)';
     v_span_index = (v_span - pv): v_span;
     
    for i=1:(m_pnts+1)
        u_span = findspan(knotU,pu,u_pnts(i));
        Nu = bspbasisDers(knotU,pu,u_pnts(i),0);
        Nu = Nu(1,:)';
       
        u_span_index = (u_span - pu): u_span;
        
        
        for j1=1:(pv+1)
            for i1=1:(pu+1)
                index = u_span_index(i1) + (v_span_index(j1)-1)*m ; 
                Z(i,j) = Z(i,j) + Uh(index)*Nu(i1)*Nv(j1);
            end
        end
        
        
    end
end

Z = Z';






z=u_Exact(X,Y);


 u_min = min(min(z));
 
 u_max = max(max(z));

 
 
% contour(X,Y,z,60,'ShowText','on');
figure(4)
pcolor(X,Y,z)
caxis([u_min,u_max])
% xlabel('x')
% ylabel('y')
set(gca,'FontName','Times New Roman','FontSize',35,'LineWidth',1);
colormap(jet)
shading interp
 dffz = max(max(abs(z-Z)))
axis equal
colorbar


figure(3)
pcolor(X,Y,Z)
caxis([u_min,u_max])
title('Contour plot of IGA solution')
set(gca,'FontName','Times New Roman','FontSize',35,'LineWidth',1);
 colormap(jet)
 shading interp
axis equal
% meshc(X,Y,Z)
% zlim([-1 0])

colorbar




        




end
