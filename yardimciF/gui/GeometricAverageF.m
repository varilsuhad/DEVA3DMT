function [Zort] = GeometricAverageF(Z)

N=length(Z);

r=abs(Z);
p=angle(Z);


% r1=power(prod(r),1/N);
% p1=sum(p)/N;


r1=geomean(r);
p1=mean(p);

Zort=r1*exp(sqrt(-1)*p1);


end

