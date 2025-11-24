% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Evaluate phase tensor variances using tipper-based inputs while guarding against invalid entries.
function [P,sigP] = variancePBF(T,e)

    al=sum(isnan(T(:)));

    if(al>0)
        P=NaN(1,2);
        sigP=P;
        return
    end

    al=sum((T(:)==0));
    if(al>0)
        P=NaN(1,2);
        sigP=P;
        return
    end

    X1=real(T(1));
    X2=real(T(2));
    Y1=imag(T(1));
    Y2=imag(T(2));

    X=[X1 X2];
    Y=[Y1 Y2];
    PSI=pinv(X)*Y;

    H=(eye(2)+PSI*PSI')^0.5;
    P=X*H;
    det=X1^2+X2^2;
    a1=X1*Y1/det;
    a2=X1*Y2/det;
    a3=X2*Y1/det;
    a4=X2*Y2/det;
    %%%%türevler  - X1
    a1tx1=(Y1*det-X1*Y1*(2*X1))/det^2;
    a2tx1=(Y2*det-X1*Y2*(2*X1))/det^2;
    a3tx1=(-X2*Y1*(2*X1))/det^2;
    a4tx1=(-X2*Y2*(2*X1))/det^2;

    %%%%türevler  - X2
    a1tx2=(-X1*Y1*(2*X2))/det^2;
    a2tx2=(-X1*Y2*(2*X2))/det^2;
    a3tx2=(Y1*det-X2*Y1*(2*X2))/det^2;
    a4tx2=(Y2*det-X2*Y2*(2*X2))/det^2;

    %%%%türevler  - Y1
    a1ty1=(X1)/det;
    a2ty1=0;
    a3ty1=X2/det;
    a4ty1=0;

    %%%%türevler  - Y2
    a1ty2=0;
    a2ty2=(X1)/det;
    a3ty2=0;
    a4ty2=X2/det;

    p1=(a1^2+a2^2+1);
    p2=(a3^2+a4^2+1);
    p3=(a1*a3+a2*a4);
    s=sqrt(p1*p2-p3^2);
    alt=sqrt(p1+p2+2*s);

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

    sigP=zeros(1,2);
    sigP(1,1)=sqrt(P1TX1^2*e(1,1)^2+P1TY1^2*e(1,1)^2+P1TX2^2*e(1,2)^2+P1TY2^2*e(1,2)^2);
    sigP(1,2)=sqrt(P2TX1^2*e(1,1)^2+P2TY1^2*e(1,1)^2+P2TX2^2*e(1,2)^2+P2TY2^2*e(1,2)^2);

end

