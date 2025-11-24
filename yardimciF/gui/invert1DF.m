% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Perform a simple 1-D inversion of apparent resistivity/phase curves.
function [ro,mf,F] = invert1DF(d,W,dz,f,set)

if(set.MT1Dfigures==1)
figure;
end
zz=cumsum(dz(:));

ro=ones(size(dz))*set.MT1Dinitalm;

lambda=set.MT1Dlambda;
lim=set.MT1Dlambdalim;
rmslim=set.MT1Drmslim;
mflim=set.MT1Dmflim;

N=length(ro);
B=[ones(N,1) -ones(N,1)*2 ones(N,1)];
A = spdiags(B,[-1 0 1],N,N);
A(1,1)=-1;A(end,end)=-1;
CC=A'*A;

W=spdiags(1./W,0,size(W,1),size(W,1));

for i=1:set.MT1Dmaxit

[Z,go,faz] = MT1DF(dz,ro,f);
if (set.MT1DInvEmp==1)
    F=Z;
else
    F=[log(go);faz];
end

dd=d-F;
mf(i)=real(dd'*W'*W*dd);
rms(i)=sqrt(mf(i)/length(dd));

J=Jacobian1DF(ro,dz,f,set);

H=J'*W'*W*J+lambda*CC;
g=J'*W'*W*dd-lambda*CC*log(ro);
dm=real(H)\real(g);

ro=exp(log(ro)+dm);

fprintf('1D inv iteration=%d misfit=%.2f rms=%.2f\n',i,mf(i),rms(i));

lambda=lambda/2;
if(lambda<=lim)
     lambda=lim;
end

if(i>2 && abs(rms(i)-rms(i-1))<set.MT1Drmschangestop)
    break
end

if(set.MT1Dfigures==1)
    plot(zz/1000,log10(ro));xlabel('z (km)');ylabel('log(ro) (\Omega.m)');
    title(strcat('iteration=',num2str(i)));
    pause(0.5);
end

if(mf(i)<=mflim || rms(i)<=rmslim)
    break;
end

end

end

