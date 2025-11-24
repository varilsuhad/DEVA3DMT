% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [PT,sigPT,PV,sigPV] = variancePTF(Z,e,T,te)

    al=sum(isnan(Z(:)));

    if(al>0)
        PT=NaN(2,2);
        sigPT=PT;
        PV=NaN(2,1);
        sigPV=PV;
        return
    end
        
    

    PT=real(Z)^-1*imag(Z);

    X=real(Z);
    detX=det(X);
    Y=imag(Z);


    sigPT=zeros(2,2);
    
    %%%P11
    a11=-X(2,2)*PT(1,1)/detX;
    a12=(-Y(2,1)+X(2,1)*PT(1,1))/detX;
    a21=X(1,2)*PT(1,1)/detX;
    a22=(Y(1,1)-X(1,1)*PT(1,1))/detX;

    b11=X(2,2)/detX;
    b21=-X(1,2)/detX;     
     
    sigPT(1,1)=sqrt(a11^2*e(1,1)^2+a21^2*e(2,1)^2+a12^2*e(1,2)^2+a22^2*e(2,2)^2 ...
                  +b11^2*e(1,1)^2+b21^2*e(2,1)^2);
    
    %%%P21    
    a11=(Y(2,1)-X(2,2)*PT(2,1))/detX;
    a12=(X(2,1)*PT(2,1))/detX;
    a21=(-Y(1,1)+X(1,2)*PT(2,1))/detX;
    a22=(-X(1,1)*PT(2,1))/detX;

    b11=-X(2,1)/detX;
    b21=X(1,1)/detX;  
    
    sigPT(2,1)=sqrt(a11^2*e(1,1)^2+a21^2*e(2,1)^2+a12^2*e(1,2)^2+a22^2*e(2,2)^2 ...
                  +b11^2*e(1,1)^2+b21^2*e(2,1)^2);
    
    %%%P12    
    a11=(-X(2,2)*PT(1,2))/detX;
    a12=(-Y(2,2)+X(2,1)*PT(1,2))/detX;
    a21=(X(1,2)*PT(1,2))/detX;
    a22=(Y(1,2)-X(1,1)*PT(1,2))/detX;

    b12=X(2,2)/detX;
    b22=-X(1,2)/detX;  
    
    sigPT(1,2)=sqrt(a11^2*e(1,1)^2+a21^2*e(2,1)^2+a12^2*e(1,2)^2+a22^2*e(2,2)^2 ...
                   +b12^2*e(1,2)^2+b22^2*e(2,2)^2);
        
    %%%P22    
    a11=(Y(2,2)-X(2,2)*PT(2,2))/detX;
    a12=(X(2,1)*PT(2,2))/detX;
    a21=(-Y(1,2)+X(1,2)*PT(2,2))/detX;
    a22=(-X(1,1)*PT(2,2))/detX;

    b12=-X(2,1)/detX;
    b22=X(1,1)/detX;  
    
    sigPT(2,2)=sqrt(a11^2*e(1,1)^2+a21^2*e(2,1)^2+a12^2*e(1,2)^2+a22^2*e(2,2)^2 ...
                  +b12^2*e(1,2)^2+b22^2*e(2,2)^2);  
              
     sigA=zeros(2,2);         
     %%%A11
     detZ=det(Z);
     a11=-Z(2,2)*Z(2,2)/detZ^2;
     a12=Z(2,2)*Z(2,1)/detZ^2;
     a21=Z(2,2)*Z(1,2)/detZ^2;
     a22=(detZ-Z(1,1)*Z(2,2))/detZ^2;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     e1=[e(1,1)^2 0;0 e(1,1)^2];
     e2=[e(1,2)^2 0;0 e(1,2)^2];
     e3=[e(2,1)^2 0;0 e(2,1)^2];
     e4=[e(2,2)^2 0;0 e(2,2)^2];
     
     q1=[a11 a12 a21 a22];
     v1=[e1 zeros(2) zeros(2) zeros(2); zeros(2) e2 zeros(2) zeros(2); ...
         zeros(2) zeros(2) e3 zeros(2); zeros(2) zeros(2) zeros(2) e4];
     
     son=q1*v1*transpose(q1);
     sigA(1,1)=sqrt(son(1,1));
     
     %%% A12
     a11=-Z(1,2)*Z(2,2)/detZ^2;
     a12=(-detZ-Z(1,2)*Z(2,1))/detZ^2;
     a21=-Z(1,2)*Z(1,2)/detZ^2;
     a22=(Z(1,2)*Z(1,1))/detZ^2;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     q1=[a11 a12 a21 a22];
     
     son=q1*v1*transpose(q1);
     sigA(1,2)=sqrt(son(1,1));    
     
     
     %%% A21
     a11=Z(2,1)*Z(2,2)/detZ^2;
     a12=(-Z(2,1)*Z(2,1))/detZ^2;
     a21=(-detZ-Z(1,2)*Z(2,1))/detZ^2;
     a22=(Z(2,1)*Z(1,1))/detZ^2;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     q1=[a11 a12 a21 a22];
     
     son=q1*v1*transpose(q1);
     sigA(2,1)=sqrt(son(1,1));      
     
     
     %%% A22
     a11=(detZ-Z(1,1)*Z(2,2))/detZ^2;
     a12=(Z(1,1)*Z(2,1))/detZ^2;
     a21=(Z(1,1)*Z(1,2))/detZ^2;
     a22=(-Z(1,1)*Z(1,1))/detZ^2;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     q1=[a11 a12 a21 a22];
     
     son=q1*v1*transpose(q1);
     sigA(2,2)=sqrt(son(1,1));     
     
     
     
     
     
     sigAT=zeros(2,1);
     te1=[te(1)^2 0;0 te(1)^2];
     te2=[te(2)^2 0;0 te(2)^2];     
     
     v1=[e1 zeros(2) zeros(2) zeros(2) zeros(2) zeros(2);...
         zeros(2) e2 zeros(2) zeros(2) zeros(2) zeros(2); ...
         zeros(2) zeros(2) e3 zeros(2) zeros(2) zeros(2); ...
         zeros(2) zeros(2) zeros(2) e4 zeros(2) zeros(2); ...
         zeros(2) zeros(2) zeros(2) zeros(2) te1 zeros(2); ...
         zeros(2) zeros(2) zeros(2) zeros(2) zeros(2) te2];
     
     %%% AT1
    
     a11=-Z(2,2)*T(1)*Z(2,2)/detZ^2 + Z(2,2)*T(2)*Z(2,1)/detZ^2;
     a12= Z(2,1)*T(1)*Z(2,2)/detZ^2 - Z(2,1)*T(2)*Z(2,1)/detZ^2;
     a21= Z(1,2)*T(1)*Z(2,2)/detZ^2 + (-T(2)*detZ-Z(1,2)*T(2)*Z(2,1))/detZ^2;
     a22=(T(1)*detZ-Z(1,1)*T(1)*Z(2,2))/detZ^2 + Z(2,1)*T(2)*Z(1,1)/detZ^2;
     b11=Z(2,2)/detZ;
     b12=-Z(2,1)/detZ;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     b11=[real(b11) -imag(b11); imag(b11) real(b11)];
     b12=[real(b12) -imag(b12); imag(b12) real(b12)];
     
     q1=[a11 a12 a21 a22 b11 b12];     
     son1=q1*v1*transpose(q1);
     sigAT(1)=sqrt(son1(2,2));
     
     %%% AT2

     a11= Z(2,2)*T(1)*Z(1,2)/detZ^2 + (T(2)*detZ-Z(1,1)*T(2)*Z(2,2))/detZ^2;
     a12=(-T(1)*detZ-Z(1,2)*T(1)*Z(2,1))/detZ^2 + Z(1,1)*T(2)*Z(2,1)/detZ^2;
     a21=-Z(1,2)*T(1)*Z(1,2)/detZ^2 + Z(1,1)*T(2)*Z(1,2)/detZ^2;
     a22= Z(1,1)*T(1)*Z(1,2)/detZ^2 - Z(1,1)*T(2)*Z(1,1)/detZ^2;     
     b11=-Z(1,2)/detZ;
     b12=Z(1,1)/detZ;
     
     a11=[real(a11) -imag(a11); imag(a11) real(a11)];
     a12=[real(a12) -imag(a12); imag(a12) real(a12)];
     a21=[real(a21) -imag(a21); imag(a21) real(a21)];
     a22=[real(a22) -imag(a22); imag(a22) real(a22)];
     b11=[real(b11) -imag(b11); imag(b11) real(b11)];
     b12=[real(b12) -imag(b12); imag(b12) real(b12)];
     
     q2=[a11 a12 a21 a22 b11 b12];     
     son2=q2*v1*transpose(q2);
     sigAT(2)=sqrt(son2(2,2));        
     
     A=inv(Z);

     
     
     K=real(A);
     G=imag(transpose(T(:))*A);
     detK=det(K);
     PV=transpose(G*inv(K));
     sigPV=zeros(2,1);
     
     %%% PT11
     a11=-K(2,2)*PV(1)/detK;
     a12= K(2,1)*PV(1)/detK;
     a21= (-G(2)+K(1,2)*PV(1))/detK;
     a22= (G(1)-K(1,1)*PV(1))/detK;
     b11= K(2,2)/detK;
     b12= -K(2,1)/detK;
     
     bak1=[a11 a12 a21 a22 b11 b12];
     
     
     sigPV(1)=sqrt(a11^2*sigA(1,1)^2+a21^2*sigA(2,1)^2+a12^2*sigA(1,2)^2+a22^2*sigA(2,2)^2 ...
                   +b11^2*sigAT(1)^2+b12^2*sigAT(2)^2);
               
     %%% PT12
     a11= (G(2)-K(2,2)*PV(2))/detK;
     a12=(-G(1)+K(2,1)*PV(2))/detK;     
     a21= K(1,2)*PV(2)/detK;
     a22=-K(1,1)*PV(2)/detK;

     b11= -K(1,2)/detK;
     b12= K(1,1)/detK;
     bak2=[a11 a12 a21 a22 b11 b12];
     
     sigPV(2)=sqrt(a11^2*sigA(1,1)^2+a21^2*sigA(2,1)^2+a12^2*sigA(1,2)^2+a22^2*sigA(2,2)^2 ...
                   +b11^2*sigAT(1)^2+b12^2*sigAT(2)^2);               

%      sigPV=sigPV/10;          
%      sigPT=sigPT/10;               
end

