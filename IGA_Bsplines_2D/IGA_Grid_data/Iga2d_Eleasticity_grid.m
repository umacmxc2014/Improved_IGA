function [Element,Coordinate,knotSpanIndex,Ubar,Vbar,m,n,Qw,NoEs,dof,uNoEs,vNoEs]=Iga2d_Eleasticity_grid(knotU,pu,knotV,pv,wights,Refinement)

% [Ubar,Vbar,Qw,dof]=RefineSurface(knotU,pu,knotV,pv,wights,Refinement);
[Qw,Ubar,Vbar,dof]=IGAknotRefineSurface(wights,knotU,pu,knotV,pv,Refinement);
[UBreaks,Uspan]=unique(Ubar); % In the old version, I use the unique to find the index of knot span;
[VBreaks,Vspan]=unique(Vbar);
% ----- Uspan(end)=Uspan(end)-pu-1;Vspan(end)=Vspan(end)-pv-1;

m=length(Ubar)-pu-1;%======= u�����ϵĻ������
n=length(Vbar)-pv-1;%======= v������ϻ������
uNoEs=length(UBreaks)-1;%==== u�����Ͻڵ������еĶϵ����
vNoEs=length(VBreaks)-1;%======== v�����ϵĽڵ������Ķϵ���



NoEs=uNoEs*vNoEs;%==== ��Ԫ����
EleDof=(pu+1)*(pv+1);%-------��Ԫ�նȾ���Ľ���
Element=zeros(NoEs,EleDof);
knotSpanIndex=zeros(NoEs,2);
Coordinate=zeros(NoEs,4);


for i1=1:uNoEs
	for j1=1:vNoEs
	row=zeros(1,EleDof);
	e=i1+(j1-1)*uNoEs;%====����������ĵ�e�ŵ�Ԫ��
	Coordinate(e,:)=[UBreaks(i1:i1+1),VBreaks(j1:j1+1)];
   % i=Uspan(i1);j=Vspan(j1);  % This is not correct for the new
   % Matlab version.
    i=findspan(Ubar,pu,UBreaks(i1)); % The knot span index of the i1-th break points in the u direction;
    j=findspan(Vbar,pv,VBreaks(j1));  % The knot span index of the j1-th break points in the v direction;
	knotSpanIndex(e,:)=[i,j];
	for k=0:pv
     temp=(k*(pu+1)+1):(k+1)*(pu+1);
	 tmp=m*(j-pv-1+k)+(i-pu:i);
     row(temp)=tmp;
	end
	Element(e,:)=row;
end
end