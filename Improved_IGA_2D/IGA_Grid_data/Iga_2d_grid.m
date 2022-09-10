function nurbsInfo=Iga_2d_grid(knotU,pu,knotV,pv,wights,Refinement)

% function [Element,Coordinate,knotSpanIndex,Ubar,Vbar,m,n,Qw,NoEs,dof]=Iga_2d_grid(knotU,pu,knotV,pv,wights,Refinement)
addpath('./NURBS')
% [Ubar,Vbar,Qw,dof]=RefineSurface(knotU,pu,knotV,pv,wights,Refinement);
[Qw,Ubar,Vbar,dof]=IGAknotRefineSurface(wights,knotU,pu,knotV,pv,Refinement);
nurbsInfo.Qw= Qw;
nurbsInfo.Ubar=Ubar;
nurbsInfo.Vbar=Vbar;
nurbsInfo.dof=dof;

UBreaks =unique(Ubar); % u 方向上节点向量中的断点的总个数.
VBreaks =unique(Vbar); % v 方向上节点向量中的断点的总个数.

nurbsInfo.UBreaks=UBreaks;
nurbsInfo.VBreaks=VBreaks;

% ----- Uspan(end)=Uspan(end)-pu-1;Vspan(end)=Vspan(end)-pv-1;

m=length(Ubar)-pu-1;%======= u 方向上基函数的个数。
n=length(Vbar)-pv-1;%======= v 方向上基函数的个数。

nurbsInfo.m=m;
nurbsInfo.n=n;
nurbsInfo.n_dofs = m*n;
nurbsInfo.pu = pu;
nurbsInfo.pv = pv;


uNoEs=length(UBreaks)-1;%==== u 方向上的区间数。
vNoEs=length(VBreaks)-1;%======== v 方向上的区间数。
NoEs=uNoEs*vNoEs;%==== 计算区域上的区间总数。

nurbsInfo.uNoEs=uNoEs;
nurbsInfo.vNoEs=vNoEs;
nurbsInfo.NoEs=NoEs;

Eledof=(pu+1)*(pv+1);%====　单元上的自由度个数。
Element=zeros(NoEs,Eledof);% 存储每个单元上的自由度的编号.
knotSpanIndex=zeros(NoEs,2);% 存储每个单元上的参数开始坐标的 knot span index. 
Coordinate=zeros(NoEs,4); % 存储每个单元上的参数坐标的四个值: $[u_i, u_{i+1}] *  [v_j, v_{j+1}]$.

for j1=1:vNoEs % 循环v方向上的全部单元.
  for i1=1:uNoEs %循环 u方向上的全部单元.
	row=zeros(1,Eledof);% 每个单元上的自由度的index.
	e=i1+(j1-1)*uNoEs;%=== 第    (i1,j1)　号单元的编号。
    % 注意，网格里单元的编号顺序是：　是先把最底部那一行单元从左到右，再从下到上来排列的.
	Coordinate(e,:)=[UBreaks(i1:i1+1),VBreaks(j1:j1+1)];% 存储当前单元的四个参数坐标.
    i=findspan(Ubar,pu,UBreaks(i1));%当前单元e上的u方向上的节点张成区间的index,即$[u_i, u_{i+1}]$.
    j=findspan(Vbar,pv,VBreaks(j1)); %当前单元e上的v方向上的节点张成区间的index,即$[v_j, v_{j+1}]$.
	knotSpanIndex(e,:)=[i,j];
	for k=0:pv
     temp=(k*(pu+1)+1):(k+1)*(pu+1);
	 tmp=m*(j-pv-1+k)+(i-pu:i);
     row(temp)=tmp;
     % 注意，这里全局自由度的编号是先把u方向上的index变化，再让v方向上的index固定，也就是:
     % N_{1,1}, N_{2,1}, ... , N_{m,1}; N_{1,2},N_{2,2},..., N_{m,2}; ...;
     % N_{1,n},N_{2,n},..., N_{m,n}.
	end
	Element(e,:)=row;
end
end




nurbsInfo.Element=Element;
nurbsInfo.Coordinate=Coordinate;
nurbsInfo.knotSpanIndex=knotSpanIndex;


bottom_row = zeros(uNoEs,pu+1);
bottom_column=zeros(uNoEs,Eledof);

top_row =bottom_row;
top_column=bottom_column;


bottom_node=zeros(uNoEs,2);
bottom_span=zeros(uNoEs,1);

for e=1:uNoEs %循环 u方向上的全部单元.  这时处理v=0或者v=1的情形.
    bottom_node(e,:)=[UBreaks(e),UBreaks(e+1)];
   
    i=findspan(Ubar,pu,UBreaks(e));%当前单元e上的u方向上的节点张成区间的index,即$[u_i, u_{i+1}]$.
    bottom_span(e)=i;
    bottom_row(e,:)=i-pu:i;
    
     j=pv+1;
     for k=0:pv
     temp=(k*(pu+1)+1):(k+1)*(pu+1);
     tmp=m*(j-pv-1+k)+(i-pu:i);
     bottom_column(e,temp)=tmp;
    end

    top_row(e,:)=m*(n-1) + (i-pu:i);  

     j=n;
     for k=0:pv
     temp=(k*(pu+1)+1):(k+1)*(pu+1);
     tmp=m*(j-pv-1+k)+(i-pu:i);
     top_column(e,temp)=tmp;
    end 
end

nurbsInfo.bottom_node=bottom_node;
nurbsInfo.bottom_span=bottom_span;

nurbsInfo.bottom_row=bottom_row;
nurbsInfo.bottom_column=bottom_column;

nurbsInfo.top_row=top_row;
nurbsInfo.top_column=top_column;




left = zeros(vNoEs,pv+1);
right=left;

for e=1:vNoEs
   j=knotSpanIndex(e,2);
   i=j-pv:j;
   left(e,:)= 1 + (i-1)*m;
   right(e,:)=m + (i-1)*m;
end 


nurbsInfo.left=left;
nurbsInfo.right=right;



bottom_dofs_2_layers =zeros(1,2*m); % The 2 layers of DOFs  near  v=0;
top_dofs_2_layers   = zeros(1,2*m);    % The 2 layers of DOFs  near  v=1;
left_dofs_2_layers   =zeros(1,2*n);
right_dofs_2_layers =zeros(1,2*n);

k=1;
    for j=1:2
        for i=1:m
        index = i + (j-1)*m;
        bottom_dofs_2_layers(k) = index;
        k = k+1;
    end
end

k=1;
    for j=(n-1):n % 对于顶端边界，我们现在先把把次外层的自由度编号排列好
        for i=1:m
        index = i + (j-1)*m;
        top_dofs_2_layers(k) = index;
        k = k+1;
    end
end




k=1;
     for j=1:n
        for i=1:2
        index = i + (j-1)*m;
        left_dofs_2_layers(k) = index;
        k = k+1;
    end
end

% left_dofs_2_layers'



k=1;
   for j=1:n 
for i=(m-1):m  % 对于右边界，我们也先把边界上的自由度编号排列好
        index = i + (j-1)*m;
        right_dofs_2_layers(k) = index;
        k = k+1;
    end
end

nurbsInfo.bottom_dofs_2_layers = bottom_dofs_2_layers;
nurbsInfo.top_dofs_2_layers = top_dofs_2_layers;
nurbsInfo.left_dofs_2_layers = left_dofs_2_layers;
nurbsInfo.right_dofs_2_layers  = right_dofs_2_layers;


% The DOF index for (i,j) is:   i + (j-1)*m 

bottom_dofs =1:m; % this is for v=0;
top_dofs   = m*(n-1) +(1:m);
left_dofs   =(0:(n-1))*m  +1;
right_dofs =(0:(n-1))*m  + m;

nurbsInfo.bottom_dofs = bottom_dofs;
nurbsInfo.top_dofs = top_dofs;

nurbsInfo.left_dofs = left_dofs;
nurbsInfo.right_dofs  = right_dofs;

bottom_dofs_2nd_layer =m + (1:m);            % v=0 边界上的往内第二层上的自由度的编号.
top_dofs_2nd_layer   = m*(n-2) +(1:m);       % v=1 边界上的往内第二层上的自由度的编号
left_dofs_2nd_layer   =(0:(n-1))*m  + 2;       % u=0 边界上的往内第二层上的自由度的编号
right_dofs_2nd_layer =(0:(n-1))*m  + m-1;   % u=1 边界上的往内第二层上的自由度的编号



nurbsInfo.bottom_dofs_2nd_layer = bottom_dofs_2nd_layer;
nurbsInfo.top_dofs_2nd_layer = top_dofs_2nd_layer;

nurbsInfo.left_dofs_2nd_layer = left_dofs_2nd_layer;
nurbsInfo.right_dofs_2nd_layer  = right_dofs_2nd_layer;
        
