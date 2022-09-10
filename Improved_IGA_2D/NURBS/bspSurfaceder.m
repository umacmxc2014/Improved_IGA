function [DSu,DSv]=bspSurfaceder(P,U,p,u,V,q,v)
%============== ����B��������(x(u,v),y(u,v),z(u,v))=S(u,v)�ڵ�(u,v)���Ĺ���u��v��һ��ƫ����ֵ��
%=========����PΪ���Ƶ���ɵ�m*n*ndim�׾���
%======== [m,n,ndim]=size(P);
%-----------  pΪu����B�����������Ĵ�����UΪ���ͽڵ�����������m+p+1���㹹�ɣ�����������Ϊm;
%%-----------  qΪv����B�����������Ĵ�����VΪ���ͽڵ�����������n+q+1���㹹�ɣ�����������Ϊn;

uspan=findspan(U,p,u);%----- u�����Ͻڵ��ų������ָ�ꣻ
vspan=findspan(V,q,v);
ndim=size(P,3);%======= ���Ƶ�����ڵĿռ�ά����Ϊ2����3��
temp=P(uspan-p:uspan,vspan-q:vspan,:);%====== �ɻ�������Ϊ0ʱ����Ӧ�Ŀ��Ƶ㹹��������u,v��һ��ƫ���Ŀ��Ƶ����
ConPtsUbar=zeros(p,q+1,ndim);%==�洢����u�����һ��ƫ��ʱ��Ҫ�Ŀ��Ƶ������ʱu����ֻҪp�����Ƶ㣻v������Ҫq+1����
ConPtsVbar=zeros(p+1,q,ndim);%==�洢����v�����һ��ƫ��ʱ��Ҫ�Ŀ��Ƶ������ʱu����ֻҪp+1�����Ƶ㣻v������Ҫq����
tempU=(U(uspan+1:uspan+p)-U(uspan-p+1:uspan))';
tempU=tempU*ones(1,q+1);%======tempUΪp*(q+1)�׾��󣬴洢u����������u��һ��ƫ��ʱ��Ҫ���¿��Ƶ�ʱ�Ķ�Ӧ�ڵ��֣�
tempV=V(vspan+1:vspan+q)-V(vspan-q+1:vspan);
tempV=ones(p+1,1)*tempV;%======= tempVΪ(p+1)*q�׾��󣬴洢v�����ϵļ������v��һ��ƫ��ʱ��Ҫ���¿��Ƶ�ڵ��֣�
for  i=1:ndim	
    ConPtsUbar(:,:,i)=p*(temp(2:end,:,i)-temp(1:end-1,:,i))./tempU;
    ConPtsVbar(:,:,i)=q*(temp(:,2:end,i)-temp(:,1:end-1,i))./tempV;
end

%% ����S(u,v)����u�����һ��ƫ������
Ubar=U;
% Ubar(1)=[];Ubar(end)=[];
 Ubar([1,end])=[];
Nu=bsplinebasis(Ubar,p-1,u);%===== ����u�������µĲ�Ϊ0��p��B���������������ڶ�u����һ�׵����������p-1�Σ�
Nv=bsplinebasis(V,q,v);
DSu=zeros(ndim,1);%====DSu�洢S(u,v)����u��һ��ƫ����ֵ��
for i=1:ndim
    temp=reshape(ConPtsUbar(:,:,i),p,q+1);
   DSu(i)=Nu'*temp*Nv;
	end

% for i=1:p
%	for j=1:q+1
%	temp=reshape(ConPtsUbar(i,j,:),ndim,1);
% DSu=DSu+temp*Nu(i)*Nv(j);
% end
% end

%% ����S(u,v)����v�����һ��ƫ������
Vbar=V;
% Vbar(1)=[];Vbar(end)=[];
Vbar([1,end])=[];
Nu=bsplinebasis(U,p,u);
Nv=bsplinebasis(Vbar,q-1,v);%===== ����v�������µĲ�Ϊ0��q��B���������������ڶ�v����һ�׵����������q-1�Σ�
DSv=zeros(ndim,1);%===== DSv�洢S(u,v)����v��һ��ƫ����ֵ��
for i=1:ndim
    temp=reshape(ConPtsVbar(:,:,i),p+1,q);
   DSv(i)=Nu'*temp*Nv;
end

% for i=1:p+1
%	for j=1:q
%	temp=reshape(ConPtsVbar(i,j,:),ndim,1);
% DSv=DSv+temp*Nu(i)*Nv(j);
% end
% end
%% end
end

%%
%=============================Test======================
% U=[0 0 0 1/2 1 1 1];p=2;u=3/10;m=4;
% V=[0 0 0 1 1 1];q=2;v=6/10;n=3;
% P=zeros(m,n,3);
% P(:,:,1)=[0 3 6 9;0 3 6 9;0 3 6 9]';
% P(:,:,2)=[0 0 0 0;2 2 2 2;4 4 4 4]';
% P(:,:,3)=[0 3 3 0;2 5 5 2;0 3 3 0]';
% [DSu,DSv]=bspSurfaceder(P,U,p,u,V,q,v)

%  ������Ϊ��
%  ������Ϊ��
% DS =
%   8.400000000000000  -0.000000000000000
%                   0   4.000000000000000
%   4.800000000000001  -0.799999999999999

%=======================================================
