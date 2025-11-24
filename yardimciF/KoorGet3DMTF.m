% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [a,b,c,err] = KoorGet3DMTF(x,y,z,xd,yd,zd)

J=zeros(3,3);

d=[xd;yd;zd];
m=[0;0;0];
err=0;
for i=1:10

[x1,y1,z1]=KoorGet3DX(x,y,z,m(1),m(2),m(3));

F=[x1;y1;z1];

dd=d-F;
mf=dd'*dd;

[J,~] = Jabc3DMTF(x,y,z,m(1),m(2),m(3),J);
J=inv(J);
J=J';

H=J'*J;
g=J'*dd;
dp=H\g;

m=dp+m;

if(mf<10^-12)
    break
end

if(i==10 && mf>10^-12)
% %     error('a,b,c kestiremedi\n');
err=-1;
end
end

a=m(1);
b=m(2);
c=m(3);

if(a<-1 || a>1 || c<-1 || c>1)
% fprintf('Yanlış var a=%f b=%f c=%f\n',a,b,c);
err=-2;
end

if(isnan(a)==1 || isnan(c)==1)
 err=-3;
end

end

