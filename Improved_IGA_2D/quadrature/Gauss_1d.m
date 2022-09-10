function s=Gauss_1d(f,np,a,b,m_row, n_column)

F=@(x)((b-a)*x + a+b)/2;
DF=(b-a)/2;

[gp,gw]=grule(np);

x=F(gp);
w=DF*gw;
s=zeros(m_row,n_column);
for i=1:np
	s=s+f(x(i))*w(i);
end

end
