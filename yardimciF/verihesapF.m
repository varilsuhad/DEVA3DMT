% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Derive MT transfer functions, apparent resistivity, and phase from field solutions at receivers.
function [veri] = verihesapF(base,x,f,nof,res,set)
mu=4*pi*10^-7;

[Ex1,Ey1,Ez1,Hx1,Hy1,Hz1] = FEhesapF(base.EL,base.NK,x(:,1),f,base.yuzey,set);  %Ey1/Hx1
[Ex2,Ey2,Ez2,Hx2,Hy2,Hz2] = FEhesapF(base.EL,base.NK,x(:,2),f,base.yuzey,set);  %Ex2/Hy2

[veri]=Zhesap3DMTF(Ex1,Ex2,Ey1,Ey2,Hx1,Hx2,Hy1,Hy2,Ez1,Ez2,Hz1,Hz2,base,nof);


ro1=abs(veri.Z).^2./(mu*2*pi*f);
tro1=abs(veri.T).^2./(mu*2*pi*f);
faz1=angle(veri.Z)/pi*180;
tfaz1=angle(veri.T)/pi*180;

veri.ro=ro1;
veri.tro=tro1;
veri.faz=faz1;
veri.tfaz=tfaz1;

veri.e1=x(:,1);
veri.e2=x(:,2);
veri.f=f;

veri.Ex1(:,1)=Ex1;
veri.Ex2(:,1)=Ex2;
veri.Ey1(:,1)=Ey1;
veri.Ey2(:,1)=Ey2;

veri.Hx1(:,1)=Hx1;
veri.Hx2(:,1)=Hx2;
veri.Hy1(:,1)=Hy1;
veri.Hy2(:,1)=Hy2;

veri.Hz1(:,1)=Hz1;
veri.Hz2(:,1)=Hz2;

veri.res=res;

if(set.saveEH==0)
    veri=rmfield(veri,'E');
    veri=rmfield(veri,'H');
end

nb=length(base.brecvlist);

% vektor=[veri.Z(:,1);...
%         veri.Z(:,2);...
%         veri.Z(:,3);...
%         veri.Z(:,4); ...
%         veri.T(:,1);...
%         veri.T(:,2);...
%         veri.PT(:,1);...
%         veri.PT(:,2);...
%         veri.PT(:,3);...
%         veri.PT(:,4);...
%         veri.PV(:,1);...
%         veri.PV(:,2);...
%         veri.PA(:,1);...
%         veri.PA(:,2);...
%         veri.PA(:,3);...
%         veri.PA(:,4);...
%         veri.PB(:,1);...
%         veri.PB(:,2);...
%         veri.V(:,1);...
%         veri.V(:,2);...
%         veri.V(:,3);...
%         veri.V(:,4);...
%         veri.V2(:,1);...
%         veri.V2(:,2);...
%         veri.V2(:,3);...
%         veri.V2(:,4)];

vektor=[veri.Z(:,1);...
        veri.Z(:,2);...
        veri.Z(:,3);...
        veri.Z(:,4); ...
        veri.T(:,1);...
        veri.T(:,2);...
        veri.PT(:,1);...
        veri.PT(:,2);...
        veri.PT(:,3);...
        veri.PT(:,4);...
        veri.PV(:,1);...
        veri.PV(:,2);...
        veri.PA(:,1);...
        veri.PA(:,2);...
        veri.PA(:,3);...
        veri.PA(:,4);...
        veri.PB(:,1);...
        veri.PB(:,2);...
        veri.Det(:,1)];

for j=1:nb
    vektor=[vektor;squeeze(veri.Qt(:,1,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Qt(:,2,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Qt(:,3,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Qt(:,4,j))];
end

for j=1:nb
    vektor=[vektor;squeeze(veri.Tt(:,1,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Tt(:,2,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Tt(:,3,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Tt(:,4,j))];
end

for j=1:nb
    vektor=[vektor;squeeze(veri.Mt(:,1,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Mt(:,2,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Mt(:,3,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Mt(:,4,j))];
end

for j=1:nb
    vektor=[vektor;squeeze(veri.Yt(:,1,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Yt(:,2,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Yt(:,3,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Yt(:,4,j))];
end

for j=1:nb
    vektor=[vektor;squeeze(veri.Ot(:,1,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Ot(:,2,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Ot(:,3,j))];
end
for j=1:nb
    vektor=[vektor;squeeze(veri.Ot(:,4,j))];
end

if(set.forwardElek==1)
    veri.Fi=base.WE.WL{nof}*vektor;
end

veri.vektor=vektor;

veri.e1=x(:,1);
veri.e2=x(:,2);

end

