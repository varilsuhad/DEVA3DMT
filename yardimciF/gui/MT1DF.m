% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [Z,roro,faz] = MT1DF(dz,ro,f)

nf=length(f);

for k=1:nf

N=length(ro);
mu=4*pi*10^-7;
w=2*pi*f(k);

for i=N:-1:1

    sig=1/ro(i);
    kn=sqrt(sqrt(-1)*w*mu*sig);

    if(i==N)
        CN=1/kn;
    else
        z=dz(i);
        ust=CN*kn+tanh(kn*z);
        alt=CN*kn*tanh(kn*z)+1;
        CN=1/kn*ust/alt;
    end

end

Z(k,1)=sqrt(-1)*w*mu*CN;

roro(k,1)=abs(Z(k))^2/(w*mu);
faz(k,1)=angle(Z(k));

end

end

