function  b=Bin(n,i)
% Revise at 22:57, 2021/11/10
% The function Bin(n,i) computes the binomial coefficient $C_{n}^i = \frac{n  (n-1) \cdots (n-i) }{i!} $
% This function is needed in DegreeElevCurve.m
numerator=1;
denominator=1;

for j=0:i-1
numerator=numerator*(n-j);
denominator=denominator*(i-j);
end
b=numerator/denominator;
