function s=Gauss_2d(f,np,a,b,c,d,m_row,n_column)
Fx=@(x)((b-a)*x + a+b)/2;
Fy=@(y)((d-c)*y+c+d)/2;
DFx=(b-a)/2;DFy=(d-c)/2;

switch np
case 1
	gp=0;gw=2;
case 2
	gp=[-0.577350269189626,0.577350269189626];gw=[1,1];
case 3
	gp=[-0.774596669241483,0,0.774596669241483];gw=[0.555555555555556,0.888888888888889,0.555555555555556];
case 4
	temp=[0.339981043584856, 0.861136311594053];tmp=[0.652145154862546,0.347854845137454];
	gp=[-temp,temp];gw=[tmp,tmp];
case 5
	temp=[0.906179845938664,0.538469310105683];tmp=[0.236926885056189,0.478628670499366];
     gp=[-temp,0,temp];gw=[tmp,0.568888888888889,tmp];
case 6
	temp=[0.932469514203152,0.661209386466265,0.238619186083197];
    tmp= [0.171324492379170,0.360761573048139,0.467913934572691];
     gp=[-temp,temp];gw=[tmp,tmp];
case 7
	temp=[0.9491079123427585,0.7415311855993945,0.4058451513773972];
    tmp=[0.1294849661688697,0.2797053914892766,0.3818300505051189];
     gp=[-temp,0,temp];gw=[tmp,0.4179591836734694,tmp];
case 8
	temp=[0.9602898564975363,0.7966664774136267,0.5255324099163290,0.1834346424956498];
    tmp=[0.1012285362903763,0.2223810344533745,0.3137066458778873,0.3626837833783620];
    gp=[-temp,temp];gw=[tmp,tmp];
  case 9
        temp=[0.9681602395076261,0.8360311073266358,0.6133714327005904,0.3242534234038089];
        gp=[-temp,0., temp];
        tmp=[0.0812743883615744,	0.1806481606948574,	0.2606106964029354,0.3123470770400029];
        gw=[tmp,0.3302393550012598,tmp];
end


s=zeros(m_row,n_column);
x=Fx(gp);wx=DFx*gw;
y=Fy(gp);wy=DFy*gw;

for i=1:np
     for j=1:np
	s=s+wx(i)*wy(j)*f(x(i),y(j));
end
end
 
end