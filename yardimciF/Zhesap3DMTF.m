% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [veri]=Zhesap3DMTF(Ex1,Ex2,Ey1,Ey2,Hx1,Hx2,Hy1,Hy2,Ez1,Ez2,Hz1,Hz2,base,ii)

if(nargin==14)

tot=length(Ex1);

for i=1:tot

       H=[Hx1(i) Hy1(i);...
          Hx2(i) Hy2(i)];
       E=[Ex1(i) Ey1(i) Hz1(i); ...
          Ex2(i) Ey2(i) Hz2(i)];

      IH=inv(H);

      O=IH*E;
      Z0(1,1)=O(1,1);  %%Zxx
      Z0(2,1)=O(1,2);  %%Zyx
      Z0(1,2)=O(2,1);  %%Zxy
      Z0(2,2)=O(2,2);  %%Zyy

      ar=base.D{ii}(i,:);
      C=[ar(1) ar(2);ar(3) ar(4)];
      CT=[ar(5) ar(6)];

      Z=C*Z0;

      T0(1,1)=O(1,3);  %%Tzx
      T0(1,2)=O(2,3);  %%Tzy

      T=CT*Z0+T0;

      PT=inv(real(Z0))*imag(Z0);
      DPT=det(real(Z0));

      A0=inv(Z0);
      PV=imag(T0*A0)*inv(real(A0));
      DPV=det(real(A0));

      IE=inv(E(1:2,1:2));

      veri.Z0(i,:)=[Z0(1,1) Z0(1,2) Z0(2,1) Z0(2,2) ];
      veri.Z(i,:)=[Z(1,1) Z(1,2) Z(2,1) Z(2,2)];
      veri.T(i,:)=[T(1) T(2)];
      veri.T0(i,:)=[T0(1) T0(2)];

      veri.PT(i,:)=[PT(1,1) PT(1,2) PT(2,1) PT(2,2) ];
      veri.PV(i,:)=[PV(1) PV(2)];

      veri.DPT(i,1)=DPT;
      veri.DPV(i,1)=DPV;
      veri.DH(i,1)=det(H);
      veri.DE(i,1)=det(E(1:2,1:2));

      veri.A0(i,:)=[A0(1,1) A0(1,2) A0(2,1) A0(2,2)];

      veri.IH(i,:)=[IH(1,1) IH(1,2) IH(2,1) IH(2,2)];
      veri.IE(i,:)=[IE(1,1) IE(1,2) IE(2,1) IE(2,2)];

      veri.TA(i,:)=T0*A0;

      veri.E(i,:)=[E(1,1) E(1,2) Ez1(i) E(2,1) E(2,2) Ez2(i)];
      veri.H(i,:)=[H(1,1) H(1,2) E(1,3) H(2,1) H(2,2) E(2,3)];

      PA0=real(Z0)*(eye(2)+PT*PT')^0.5;
      veri.PA0(i,:)=[PA0(1,1) PA0(1,2) PA0(2,1) PA0(2,2)];
      PA=C*PA0;
      veri.PA(i,:)=[PA(1,1) PA(1,2) PA(2,1) PA(2,2)];

      PTV=pinv(real(T))*imag(T);
      PB=real(T)*(eye(2)+PTV*PTV')^0.5;
      veri.PTV(i,:)=[PTV(1,1) PTV(1,2) PTV(2,1) PTV(2,2)];

      veri.PB(i,:)=PB;
      veri.DPB(i,1)=real(T)*real(T)';

      veri.Det(i,1)=sqrt((Z(1,1)^2+Z(1,2)^2+Z(2,1)^2+Z(2,2)^2)/2);





end

lst=base.brecvlist;
veri.IEb=zeros(0,4);
veri.IHb=zeros(0,4);
nb=length(lst);

veri.Tt0=zeros(0,4,nb);
veri.Qt0=zeros(0,4,nb);
veri.Tt=zeros(0,4,nb);
veri.Qt=zeros(0,4,nb);
veri.Mt=zeros(0,4,nb);
veri.Ot=zeros(0,4,nb);
veri.Yt=zeros(0,4,nb);

for j=1:length(lst)
    jj=lst(j);
    Hb=[Hx1(jj) Hx2(jj);...
       Hy1(jj) Hy2(jj)];
    Eb=[Ex1(jj) Ex2(jj); ...
       Ey1(jj) Ey2(jj)];

    IEb=inv(Eb);
    IHb=inv(Hb);

    veri.IEb(j,:)=[IEb(1,1) IEb(2,1) IEb(1,2) IEb(2,2)];
    veri.IHb(j,:)=[IHb(1,1) IHb(2,1) IHb(1,2) IHb(2,2)];

    for i=1:tot
    H=[Hx1(i) Hx2(i);...
       Hy1(i) Hy2(i)];
    E=[Ex1(i) Ex2(i); ...
       Ey1(i) Ey2(i)];

    ar=base.D{ii}(i,:);
    C=[ar(1) ar(2);ar(3) ar(4)];

    Q0=E*(IHb);
    T0=E*(IEb);
    M=H*(IHb);

    Q=C*Q0;
    T=C*T0;

    Y=inv(real(Q))*imag(Q);
    O=inv(real(T))*imag(T);

    veri.detQ(j,i)=det(real(Q0));
    veri.detT(j,i)=det(real(T0));

    veri.Qt0(i,:,j)=[Q0(1,1) Q0(1,2) Q0(2,1) Q0(2,2)];
    veri.Tt0(i,:,j)=[T0(1,1) T0(1,2) T0(2,1) T0(2,2)];

    veri.Qt(i,:,j)=[Q(1,1) Q(1,2) Q(2,1) Q(2,2)];
    veri.Tt(i,:,j)=[T(1,1) T(1,2) T(2,1) T(2,2)];
    veri.Mt(i,:,j)=[M(1,1) M(1,2) M(2,1) M(2,2)];

    veri.Yt(i,:,j)=[Y(1,1) Y(1,2) Y(2,1) Y(2,2)];
    veri.Ot(i,:,j)=[O(1,1) O(1,2) O(2,1) O(2,2)];
    end
end

else

tot=length(Ex1);

for i=1:tot

       H=[Hx1(i) Hy1(i);...
          Hx2(i) Hy2(i)];
       E=[Ex1(i) Ey1(i) Hz1(i); ...
          Ex2(i) Ey2(i) Hz2(i)];

      IH=inv(H);

      O=IH*E;
      Z0(1,1)=O(1,1);  %%Zxx
      Z0(2,1)=O(1,2);  %%Zyx
      Z0(1,2)=O(2,1);  %%Zxy
      Z0(2,2)=O(2,2);  %%Zyy
      T0(1,1)=O(1,3);  %%Tzx
      T0(1,2)=O(2,3);  %%Tzy

      veri(i,:)=[Z0(1,1) Z0(2,1) Z0(1,2) Z0(2,2) T0(1,1) T0(1,2)];

end

end

end
