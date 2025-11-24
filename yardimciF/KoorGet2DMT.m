% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Map a physical quadrilateral (x,z) location to local coordinates within the reference element.
function [a,c,err] = KoorGet2DMT(x,z,xd,zd)

J=zeros(2,2);

d=[xd;zd];
m=[0;0];
err=0;
for i=1:10

% [x1,y1,z1]=KoorGetX(x,y,z,m(1),m(2),m(3));
[x1,z1] = KoorGet2DX(x,z,m(1),m(2));

F=[x1;z1];

dd=d-F;
mf=dd'*dd;

[J,det1] = Jabc2DMT(x,z,m(1),m(2),J,1);
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
c=m(2);

if(a<-1 || a>1 || c<-1 || c>1)
% fprintf('Yanlış var a=%f b=%f c=%f\n',a,b,c);
err=-2;
end

if(isnan(a)==1 || isnan(c)==1)
 err=-3;
end

end

