function plotGrid_Surface(ConPts,weights,knotU,pu,knotV,pv,Refinement)


 [~,Ubar,Vbar,~]=IGAknotRefineSurface(weights,knotU,pu,knotV,pv,Refinement);

Resolution=0.01/4;
u=knotU(1):Resolution:knotU(end);nu=length(u);
v=knotV(1):Resolution:knotV(end);nv=length(v);
ndim=size(ConPts,3);
UBreaks=unique(Ubar);NoUBreaks=length(UBreaks);
VBreaks=unique(Vbar);NoVBreaks=length(VBreaks);
SBreaks=zeros(NoUBreaks,nv,ndim);



 for i=1:NoUBreaks
     for j=1:nv
     SBreaks(i,j,:)=PointOnNurbsSurface(ConPts,weights,knotU,pu,UBreaks(i),knotV,pv,v(j));
  end
      x=reshape(SBreaks(i,:,1),nv,1);
      y=reshape(SBreaks(i,:,2),nv,1);

      plot(x,y,'k','LineWidth',3)
            
    
      
 hold on
 axis equal
set(gca,'FontName','Times New Roman','FontSize',35,'LineWidth',2);

end
hold on
SBreaks=zeros(nu,NoVBreaks,ndim);

	 for j=1:NoVBreaks
     for i=1:nu
	  SBreaks(i,j,:)=PointOnNurbsSurface(ConPts,weights,knotU,pu,u(i),knotV,pv,VBreaks(j));
  end
      x=reshape(SBreaks(:,j,1),nu,1);
      y=reshape(SBreaks(:,j,2),nu,1);
      
      plot(x,y,'k','LineWidth',3)

      hold on
     end 

axis equal
 hold off 

