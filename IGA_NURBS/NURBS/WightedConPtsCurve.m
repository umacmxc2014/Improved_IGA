function Pw=WightedConPtsCurve(P,w)
%===�����NUBRS��������������ʾʱ����Ҫ�õ��Ĵ�Ȩϵ���Ŀ��Ƶ���󣻿��Ƶ����PΪndim*n��;
%== nΪ���Ƶ�ĸ�������u�����ϻ������ĸ���Ϊn��
%=== �����ΪPw��һ��(ndim+1)*n�͵Ĵ�Ȩϵ���Ŀ��Ƶ����

[ndim,n]=size(P);
Pw=zeros(ndim+1,n);
Pw(end,:)=w;
for i=1:ndim
	Pw(i,:)=P(i,:).*w;
end
    
end