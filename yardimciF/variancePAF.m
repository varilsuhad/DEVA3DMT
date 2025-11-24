% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [P,sigP] = variancePAF(Z,e)

    al=sum(isnan(Z(:)));

    if(al>0)
        P=NaN(2,2);
        sigP=P;
        return
    end

    X=real(Z);
    Y=imag(Z);
    PHI=inv(X)*Y;

    X1=X(1,1);
    X2=X(1,2);
    X3=X(2,1);
    X4=X(2,2);

    Y1=Y(1,1);
    Y2=Y(1,2);
    Y3=Y(2,1);
    Y4=Y(2,2);

    det=X1*X4-X2*X3;
    a1=(X4*Y1-X2*Y3)/det;
    a2=(X4*Y2-X2*Y4)/det;
    a3=(X1*Y3-X3*Y1)/det;
    a4=(X1*Y4-X3*Y2)/det;

%%%%türevler  - X1
a1tx1=-X4*(X4*Y1-X2*Y3)/det^2;
a2tx1=-X4*(X4*Y2-X2*Y4)/det^2;
a3tx1=(Y3*det-X4*(X1*Y3-X3*Y1))/det^2;
a4tx1=(Y4*det-X4*(X1*Y4-X3*Y2))/det^2;
%%%%türevler - X2
a1tx2=((X4*Y1-X2*Y3)*X3-Y3*det)/det^2;
a2tx2=((X4*Y2-X2*Y4)*X3-Y4*det)/det^2;
a3tx2=(X1*Y3-X3*Y1)*X3/det^2;
a4tx2=(X1*Y4-X3*Y2)*X3/det^2;
%%%%türevler  - X3
a1tx3=(X4*Y1-X2*Y3)*X2/det^2;
a2tx3=(X4*Y2-X2*Y4)*X2/det^2;
a3tx3=(X2*(X1*Y3-X3*Y1)-Y1*det)/det^2;
a4tx3=(X2*(X1*Y4-X3*Y2)-Y2*det)/det^2;
%%%%türevler - X4
a1tx4=(Y1*det-(X4*Y1-X2*Y3)*X1)/det^2;
a2tx4=(Y2*det-(X4*Y2-X2*Y4)*X1)/det^2;
a3tx4=-X1*(X1*Y3-X3*Y1)/det^2;
a4tx4=-X1*(X1*Y4-X3*Y2)/det^2;
%%%%türevler - Y1
a1ty1=X4/det;
a2ty1=0;
a3ty1=-X3/det;
a4ty1=0;
%%%%türevler  - Y2
a1ty2=0;
a2ty2=X4/det;
a3ty2=0;
a4ty2=-X3/det;
%%%%türevler- Y3
a1ty3=-X2/det;
a2ty3=0;
a3ty3=X1/det;
a4ty3=0;
%Türevler Y4'e göre
a1ty4=0;
a2ty4=-X2/det;
a3ty4=0;
a4ty4=X1/det;

%bunlar sabit
p1=(a1^2+a2^2+1);
p2=(a3^2+a4^2+1);
p3=(a1*a3+a2*a4);
s=sqrt(p1*p2-p3^2);

alt=sqrt(p1+p2+2*s);

P=zeros(2,2);

P(1,1)=(X1*(p1+s)+X2*p3)/alt;
P(1,2)=(X1*p3+X2*(p2+s))/alt;
P(2,1)=(X3*(p1+s)+X4*p3)/alt;
P(2,2)=(X3*p3+X4*(p2+s))/alt;

%P1 - X1
p1t=2*a1*a1tx1+2*a2*a2tx1;
p2t=2*a3*a3tx1+2*a4*a4tx1;
p3t=a1*a3tx1+a1tx1*a3+a2*a4tx1+a4*a2tx1;

ust=X1*(p1+s)+X2*p3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(p1+s)+X1*(p1t+stur)+X2*p3t;
P1TX1=(alt*ustt-altt*ust)/alt^2;

%P1 - X2
p1t=2*a1*a1tx2+2*a2*a2tx2;
p2t=2*a3*a3tx2+2*a4*a4tx2;
p3t=a1*a3tx2+a1tx2*a3+a2*a4tx2+a4*a2tx2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t+p3;
P1TX2=(alt*ustt-altt*ust)/alt^2;

%P1 - X3
p1t=2*a1*a1tx3+2*a2*a2tx3;
p2t=2*a3*a3tx3+2*a4*a4tx3;
p3t=a1*a3tx3+a1tx3*a3+a2*a4tx3+a4*a2tx3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TX3=(alt*ustt-altt*ust)/alt^2;

%P1 - X4
p1t=2*a1*a1tx4+2*a2*a2tx4;
p2t=2*a3*a3tx4+2*a4*a4tx4;
p3t=a1*a3tx4+a1tx4*a3+a2*a4tx4+a4*a2tx4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TX4=(alt*ustt-altt*ust)/alt^2;

%P1 - Y1
p1t=2*a1*a1ty1+2*a2*a2ty1;
p2t=2*a3*a3ty1+2*a4*a4ty1;
p3t=a1*a3ty1+a1ty1*a3+a2*a4ty1+a4*a2ty1;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TY1=(alt*ustt-altt*ust)/alt^2;

%P1 - Y2
p1t=2*a1*a1ty2+2*a2*a2ty2;
p2t=2*a3*a3ty2+2*a4*a4ty2;
p3t=a1*a3ty2+a1ty2*a3+a2*a4ty2+a4*a2ty2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TY2=(alt*ustt-altt*ust)/alt^2;

%P1 - Y3
p1t=2*a1*a1ty3+2*a2*a2ty3;
p2t=2*a3*a3ty3+2*a4*a4ty3;
p3t=a1*a3ty3+a1ty3*a3+a2*a4ty3+a4*a2ty3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TY3=(alt*ustt-altt*ust)/alt^2;

%P1 - Y4
p1t=2*a1*a1ty4+2*a2*a2ty4;
p2t=2*a3*a3ty4+2*a4*a4ty4;
p3t=a1*a3ty4+a1ty4*a3+a2*a4ty4+a4*a2ty4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*(p1t+stur)+X2*p3t;
P1TY4=(alt*ustt-altt*ust)/alt^2;

%P2 - X1
p1t=2*a1*a1tx1+2*a2*a2tx1;
p2t=2*a3*a3tx1+2*a4*a4tx1;
p3t=a1*a3tx1+a1tx1*a3+a2*a4tx1+a4*a2tx1;

ust=(X1*p3+X2*(p2+s));

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=p3+X1*p3t+X2*(p2t+stur);
P2TX1=(alt*ustt-altt*ust)/alt^2;

%P2 - X2
p1t=2*a1*a1tx2+2*a2*a2tx2;
p2t=2*a3*a3tx2+2*a4*a4tx2;
p3t=a1*a3tx2+a1tx2*a3+a2*a4tx2+a4*a2tx2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur)+(p2+s);
P2TX2=(alt*ustt-altt*ust)/alt^2;

%P2 - X3
p1t=2*a1*a1tx3+2*a2*a2tx3;
p2t=2*a3*a3tx3+2*a4*a4tx3;
p3t=a1*a3tx3+a1tx3*a3+a2*a4tx3+a4*a2tx3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TX3=(alt*ustt-altt*ust)/alt^2;

%P2 - X4
p1t=2*a1*a1tx4+2*a2*a2tx4;
p2t=2*a3*a3tx4+2*a4*a4tx4;
p3t=a1*a3tx4+a1tx4*a3+a2*a4tx4+a4*a2tx4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TX4=(alt*ustt-altt*ust)/alt^2;

%P2 - Y1
p1t=2*a1*a1ty1+2*a2*a2ty1;
p2t=2*a3*a3ty1+2*a4*a4ty1;
p3t=a1*a3ty1+a1ty1*a3+a2*a4ty1+a4*a2ty1;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TY1=(alt*ustt-altt*ust)/alt^2;

%P2 - Y2
p1t=2*a1*a1ty2+2*a2*a2ty2;
p2t=2*a3*a3ty2+2*a4*a4ty2;
p3t=a1*a3ty2+a1ty2*a3+a2*a4ty2+a4*a2ty2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TY2=(alt*ustt-altt*ust)/alt^2;

%P2 - Y3
p1t=2*a1*a1ty3+2*a2*a2ty3;
p2t=2*a3*a3ty3+2*a4*a4ty3;
p3t=a1*a3ty3+a1ty3*a3+a2*a4ty3+a4*a2ty3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TY3=(alt*ustt-altt*ust)/alt^2;

%P2 - Y4
p1t=2*a1*a1ty4+2*a2*a2ty4;
p2t=2*a3*a3ty4+2*a4*a4ty4;
p3t=a1*a3ty4+a1ty4*a3+a2*a4ty4+a4*a2ty4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X1*p3t+X2*(p2t+stur);
P2TY4=(alt*ustt-altt*ust)/alt^2;

%P3 - X1
p1t=2*a1*a1tx1+2*a2*a2tx1;
p2t=2*a3*a3tx1+2*a4*a4tx1;
p3t=a1*a3tx1+a1tx1*a3+a2*a4tx1+a4*a2tx1;

ust=(X3*(p1+s)+X4*p3);

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TX1=(alt*ustt-altt*ust)/alt^2;

%P3 - X2
p1t=2*a1*a1tx2+2*a2*a2tx2;
p2t=2*a3*a3tx2+2*a4*a4tx2;
p3t=a1*a3tx2+a1tx2*a3+a2*a4tx2+a4*a2tx2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TX2=(alt*ustt-altt*ust)/alt^2;

%P3 - X3
p1t=2*a1*a1tx3+2*a2*a2tx3;
p2t=2*a3*a3tx3+2*a4*a4tx3;
p3t=a1*a3tx3+a1tx3*a3+a2*a4tx3+a4*a2tx3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+(p1+s)+X4*p3t);
P3TX3=(alt*ustt-altt*ust)/alt^2;

%P3 - X4
p1t=2*a1*a1tx4+2*a2*a2tx4;
p2t=2*a3*a3tx4+2*a4*a4tx4;
p3t=a1*a3tx4+a1tx4*a3+a2*a4tx4+a4*a2tx4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t+p3);
P3TX4=(alt*ustt-altt*ust)/alt^2;

%P3 - Y1
p1t=2*a1*a1ty1+2*a2*a2ty1;
p2t=2*a3*a3ty1+2*a4*a4ty1;
p3t=a1*a3ty1+a1ty1*a3+a2*a4ty1+a4*a2ty1;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TY1=(alt*ustt-altt*ust)/alt^2;

%P3 - Y2
p1t=2*a1*a1ty2+2*a2*a2ty2;
p2t=2*a3*a3ty2+2*a4*a4ty2;
p3t=a1*a3ty2+a1ty2*a3+a2*a4ty2+a4*a2ty2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TY2=(alt*ustt-altt*ust)/alt^2;

%P3 - Y3
p1t=2*a1*a1ty3+2*a2*a2ty3;
p2t=2*a3*a3ty3+2*a4*a4ty3;
p3t=a1*a3ty3+a1ty3*a3+a2*a4ty3+a4*a2ty3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TY3=(alt*ustt-altt*ust)/alt^2;

%P3 - Y4
p1t=2*a1*a1ty4+2*a2*a2ty4;
p2t=2*a3*a3ty4+2*a4*a4ty4;
p3t=a1*a3ty4+a1ty4*a3+a2*a4ty4+a4*a2ty4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=(X3*(p1t+stur)+X4*p3t);
P3TY4=(alt*ustt-altt*ust)/alt^2;

%P4 - X1
p1t=2*a1*a1tx1+2*a2*a2tx1;
p2t=2*a3*a3tx1+2*a4*a4tx1;
p3t=a1*a3tx1+a1tx1*a3+a2*a4tx1+a4*a2tx1;

ust=(X3*p3+X4*(p2+s));

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TX1=(alt*ustt-altt*ust)/alt^2;

%P4 - X2
p1t=2*a1*a1tx2+2*a2*a2tx2;
p2t=2*a3*a3tx2+2*a4*a4tx2;
p3t=a1*a3tx2+a1tx2*a3+a2*a4tx2+a4*a2tx2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TX2=(alt*ustt-altt*ust)/alt^2;

%P4 - X3
p1t=2*a1*a1tx3+2*a2*a2tx3;
p2t=2*a3*a3tx3+2*a4*a4tx3;
p3t=a1*a3tx3+a1tx3*a3+a2*a4tx3+a4*a2tx3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+p3+X4*(p2t+stur);
P4TX3=(alt*ustt-altt*ust)/alt^2;

%P4 - X4
p1t=2*a1*a1tx4+2*a2*a2tx4;
p2t=2*a3*a3tx4+2*a4*a4tx4;
p3t=a1*a3tx4+a1tx4*a3+a2*a4tx4+a4*a2tx4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur)+(p2+s);
P4TX4=(alt*ustt-altt*ust)/alt^2;

%P4 - Y1
p1t=2*a1*a1ty1+2*a2*a2ty1;
p2t=2*a3*a3ty1+2*a4*a4ty1;
p3t=a1*a3ty1+a1ty1*a3+a2*a4ty1+a4*a2ty1;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TY1=(alt*ustt-altt*ust)/alt^2;

%P4 - Y2
p1t=2*a1*a1ty2+2*a2*a2ty2;
p2t=2*a3*a3ty2+2*a4*a4ty2;
p3t=a1*a3ty2+a1ty2*a3+a2*a4ty2+a4*a2ty2;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TY2=(alt*ustt-altt*ust)/alt^2;

%P4 - Y3
p1t=2*a1*a1ty3+2*a2*a2ty3;
p2t=2*a3*a3ty3+2*a4*a4ty3;
p3t=a1*a3ty3+a1ty3*a3+a2*a4ty3+a4*a2ty3;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TY3=(alt*ustt-altt*ust)/alt^2;

%P4 - Y4
p1t=2*a1*a1ty4+2*a2*a2ty4;
p2t=2*a3*a3ty4+2*a4*a4ty4;
p3t=a1*a3ty4+a1ty4*a3+a2*a4ty4+a4*a2ty4;

stur=0.5*inv(s)*(p1t*p2+p1*p2t-2*p3*p3t);
altt=0.5*inv(alt)*(p1t+p2t+2*stur);
ustt=X3*p3t+X4*(p2t+stur);
P4TY4=(alt*ustt-altt*ust)/alt^2;

sigP=zeros(2,2);

sigP(1,1)=sqrt(P1TX1^2*e(1,1)^2+P1TY1^2*e(1,1)^2+P1TX2^2*e(1,2)^2+P1TY2^2*e(1,2)^2 ...
                +P1TX3^2*e(2,1)^2+P1TY3^2*e(2,1)^2+P1TX4^2*e(2,1)^2+P1TY4^2*e(2,1)^2);

sigP(1,2)=sqrt(P2TX1^2*e(1,1)^2+P2TY1^2*e(1,1)^2+P2TX2^2*e(1,2)^2+P2TY2^2*e(1,2)^2 ...
                +P2TX3^2*e(2,1)^2+P2TY3^2*e(2,1)^2+P2TX4^2*e(2,2)^2+P2TY4^2*e(2,2)^2);

sigP(2,1)=sqrt(P3TX1^2*e(1,1)^2+P3TY1^2*e(1,1)^2+P3TX2^2*e(1,2)^2+P3TY2^2*e(1,2)^2 ...
                +P3TX3^2*e(2,1)^2+P3TY3^2*e(2,1)^2+P3TX4^2*e(2,2)^2+P3TY4^2*e(2,2)^2);

sigP(2,2)=sqrt(P4TX1^2*e(1,1)^2+P4TY1^2*e(1,1)^2+P4TX2^2*e(1,2)^2+P4TY2^2*e(1,2)^2 ...
                +P4TX3^2*e(2,1)^2+P4TY3^2*e(2,1)^2+P4TX4^2*e(2,2)^2+P4TY4^2*e(2,2)^2);

end

