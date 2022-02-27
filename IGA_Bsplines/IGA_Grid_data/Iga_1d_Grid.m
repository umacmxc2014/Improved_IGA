function [Element,Coordinate,Ubar,wbar,NoEs,dof]=Iga_1d_Grid(U,w,p,Refinement)

 [Ubar,wbar]=IGAknotRefineCurve(U,w,p,Refinement);
 m=length(Ubar);
 n=m-p-1;
 dof=n;
 [UBreaks]=unique(Ubar);
 NoBreaks=length(UBreaks);
 NoEs=NoBreaks-1;
 dof_e=p+1;
 Element=zeros(NoEs,dof_e);Coordinate=zeros(NoEs,2);
 for e=1:NoEs
  i=findspan(Ubar,p,UBreaks(e));% In the last version code, I use the  unique function to find knot span; 
 Element(e,:)=i-p:i;
 Coordinate(e,:)=UBreaks(e:e+1);
 end
 