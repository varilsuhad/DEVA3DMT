% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Filter unusable sites/components and assemble weighting for all data classes before inversion.
function [WE,totD,totC,base] = elimination3DMTALLFFF(data,data2,set,base)

[nf,nist,~]=size(data);
[nf,nist,~,nb]=size(data2);

katWZ=set.katWZ;
katWW=set.katWW;
katWPT=set.katWPT;
katWPV=set.katWPV;
katWQ=set.katWQ;
katWT=set.katWT;
katWM=set.katWM;
katWY=set.katWY;
katWO=set.katWO;
katPA=set.katPA;
katPB=set.katPB;
katWD=set.katWD;

tres=10^-5;
c2=0;

%%%Data ele
for i=1:nf
    c=0;
    cc=0;c3=0;
    clear ix iy iv ww cp
    ix=[];iy=[];iv=[];ww=[];

    %%% Zxx
    for k=1:nist
    tk1=data(i,k,1);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(1)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,7)*katWZ;
    c2=c2+1;wwD(c2,1)=1./data(i,k,7)*katWZ;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% Zxy
    for k=1:nist
    tk1=data(i,k,2);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(2)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,8)*katWZ;
    c2=c2+1;wwD(c2,1)=1./data(i,k,8)*katWZ;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% Zyx
    for k=1:nist
    tk1=data(i,k,3);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(3)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,9)*katWZ;
    c2=c2+1;wwD(c2,1)=1./data(i,k,9)*katWZ;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% Zyy
    for k=1:nist
    tk1=data(i,k,4);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(4)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,10)*katWZ;
    c2=c2+1;wwD(c2,1)=1./data(i,k,10)*katWZ;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% Tzx
    for k=1:nist
    tk1=data(i,k,5);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(7)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,11)*katWW;
    c2=c2+1;wwD(c2,1)=1./data(i,k,11)*katWW;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% Tzy
    for k=1:nist
    tk1=data(i,k,6);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw(8)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data(i,k,12)*katWW;
    c2=c2+1;wwD(c2,1)=1./data(i,k,12)*katWW;
    cp(c,1)=ww(c,1)*tk1;
    end

    %%% PTxx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,13)=0;
    data(i,k,19)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(11)=='0')
    continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);
    if(set.rationalvar==1)
    sigPT=ones(size(sigPT))*max(PT(:))*0.05;
    sigPV=ones(size(sigPV))*max(PV(:))*0.05;
    end

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,1)*katWPT;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,1)*katWPT;
    cp(c,1)=ww(c,1)*PT(1,1);
    data(i,k,13)=PT(1,1);
    data(i,k,19)=sigPT(1,1);
    end

    %%% PTxy
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,14)=0;
    data(i,k,20)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(12)=='0')
    continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);
    if(set.rationalvar==1)
    sigPT=ones(size(sigPT))*max(PT(:))*0.05;
    sigPV=ones(size(sigPV))*max(PV(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,2)*katWPT;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,2)*katWPT;
    cp(c,1)=ww(c,1)*PT(1,2);
    data(i,k,14)=PT(1,2);
    data(i,k,20)=sigPT(1,2);
    end

    %%% PTyx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,15)=0;
    data(i,k,21)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(13)=='0')
    continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);
    if(set.rationalvar==1)
    sigPT=ones(size(sigPT))*max(PT(:))*0.05;
    sigPV=ones(size(sigPV))*max(PV(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,1)*katWPT;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,1)*katWPT;
    cp(c,1)=ww(c,1)*PT(2,1);
    data(i,k,15)=PT(2,1);
    data(i,k,21)=sigPT(2,1);
    end

    %%% PTyy
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,16)=0;
    data(i,k,22)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(14)=='0')
        continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);
    if(set.rationalvar==1)
    sigPT=ones(size(sigPT))*max(PT(:))*0.05;
    sigPV=ones(size(sigPV))*max(PV(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,2)*katWPT;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,2)*katWPT;
    cp(c,1)=ww(c,1)*PT(2,2);
    data(i,k,16)=PT(2,2);
    data(i,k,22)=sigPT(2,2);
    end

    %%% PVzx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,17)=0;
    data(i,k,23)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || ...
           isnan(tk5)==1 || isnan(tk6)==1 || set.sw(15)=='0')
        continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPV(1)*katWPV;
    c2=c2+1;wwD(c2,1)=1./sigPV(1)*katWPV;
    cp(c,1)=ww(c,1)*PV(1);
    data(i,k,17)=PV(1);
    data(i,k,23)=sigPV(1);
    end

    %%% PVzy
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);
    tk5=data(i,k,5);
    tk6=data(i,k,6);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);
    te5=data(i,k,11);
    te6=data(i,k,12);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    T=[tk5;tk6];
    Te=[te5;te6];
    data(i,k,18)=0;
    data(i,k,24)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || ...
           isnan(tk5)==1 || isnan(tk6)==1 || set.sw(16)=='0')
        continue;
    end
    [PT,sigPT,PV,sigPV] = variancePTF(Z,Ze,T,Te);
    if(set.rationalvar==1)
    sigPT=ones(size(sigPT))*max(PT(:))*0.05;
    sigPV=ones(size(sigPV))*max(PV(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPV(2)*katWPV;
    c2=c2+1;wwD(c2,1)=1./sigPV(2)*katWPV;
    cp(c,1)=ww(c,1)*PV(2);
    data(i,k,18)=PV(2);
    data(i,k,24)=sigPV(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%AMPLITUDE TENSOR%%%%%%%%%%%%%%%%%%%%%%
    %%% Pxx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data(i,k,25)=0;
    data(i,k,29)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(18)=='0')
        continue;
    end
    [PA,sigPA] = variancePAF(Z,Ze);
    if(set.rationalvar==1)
    sigPA=ones(size(sigPA))*max(PA(:))*0.05;
    end

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPA(1,1)*katPA;
    c2=c2+1;wwD(c2,1)=1./sigPA(1,1)*katPA;
    cp(c,1)=ww(c,1)*PA(1,1);

    data(i,k,25)=PA(1,1);
    data(i,k,29)=sigPA(1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%% Pxy
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data(i,k,26)=0;
    data(i,k,30)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(19)=='0')
        continue;
    end
    [PA,sigPA] = variancePAF(Z,Ze);
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPA(1,2)*katPA;
    c2=c2+1;wwD(c2,1)=1./sigPA(1,2)*katPA;
    cp(c,1)=ww(c,1)*PA(1,2);

    data(i,k,26)=PA(1,2);
    data(i,k,30)=sigPA(1,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%% Pyx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data(i,k,27)=0;
    data(i,k,31)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(20)=='0')
        continue;
    end
    [PA,sigPA] = variancePAF(Z,Ze);
    if(set.rationalvar==1)
    sigPA=ones(size(sigPA))*max(PA(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPA(2,1)*katPA;
    c2=c2+1;wwD(c2,1)=1./sigPA(2,1)*katPA;
    cp(c,1)=ww(c,1)*PA(2,1);

    data(i,k,27)=PA(2,1);
    data(i,k,31)=sigPA(2,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%% Pyy
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data(i,k,28)=0;
    data(i,k,32)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(21)=='0')
        continue;
    end
    [PA,sigPA] = variancePAF(Z,Ze);
    if(set.rationalvar==1)
    sigPA=ones(size(sigPA))*max(PA(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPA(2,2)*katPA;
    c2=c2+1;wwD(c2,1)=1./sigPA(2,2)*katPA;
    cp(c,1)=ww(c,1)*PA(2,2);

    data(i,k,28)=PA(2,2);
    data(i,k,32)=sigPA(2,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%AMPLITUDE VECTOR%%%%%%%%%%%%%%%%%%%%%%

    %%% Pzx
    for k=1:nist
    tk1=data(i,k,5);
    tk2=data(i,k,6);

    te1=data(i,k,11);
    te2=data(i,k,12);

    T=[tk1 tk2];
    Te=[te1 te2];

    data(i,k,33)=0;
    data(i,k,35)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || (tk1)==0 || (tk2)==0 || set.sw(22)=='0')
        continue;
    end
    [PB,sigPB] = variancePBF(T,Te);
    if(set.rationalvar==1)
    sigPB=ones(size(sigPB))*max(PB(:))*0.05;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPB(1)*katPB;
    c2=c2+1;wwD(c2,1)=1./sigPB(1)*katPB;
    cp(c,1)=ww(c,1)*PB(1);

    data(i,k,33)=PB(1);
    data(i,k,35)=sigPB(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%% Pzy
    for k=1:nist
    tk1=data(i,k,5);
    tk2=data(i,k,6);

    te1=data(i,k,11);
    te2=data(i,k,12);

    T=[tk1 tk2];
    Te=[te1 te2];

    data(i,k,34)=0;
    data(i,k,36)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || (tk1)==0 || (tk2)==0 || set.sw(23)=='0')
        continue;
    end
    [PB,sigPB] = variancePBF(T,Te);
    if(set.rationalvar==1)
    sigPB=ones(size(sigPB))*max(PB(:))*0.05;
    end

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPB(2)*katPB;
    c2=c2+1;wwD(c2,1)=1./sigPB(2)*katPB;
    cp(c,1)=ww(c,1)*PB(2);

    data(i,k,34)=PB(2);
    data(i,k,36)=sigPB(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%Determinant%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dxx
    for k=1:nist
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    tk3=data(i,k,3);
    tk4=data(i,k,4);

    te1=data(i,k,7);
    te2=data(i,k,8);
    te3=data(i,k,9);
    te4=data(i,k,10);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data(i,k,37)=0;
    data(i,k,38)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw(25)=='0')
        continue;
    end
    [Det,sigD] = varianceDF(Z,Ze);
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigD(1,1)*katWD;
    c2=c2+1;wwD(c2,1)=1./sigD(1,1)*katWD;
    cp(c,1)=ww(c,1)*Det(1,1);

    data(i,k,37)=Det(1,1);
    data(i,k,38)=sigD(1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%INTER-SITE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Qxx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,1,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(1)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,13,j)*katWQ;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,13,j)*katWQ;
    end
    end

    %%% Qxy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,2,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(2)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,14,j)*katWQ;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,14,j)*katWQ;
    end
    end

    %%% Qyx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,3,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(3)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,15,j)*katWQ;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,15,j)*katWQ;
    end
    end

    %%% Qyy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,4,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(4)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,16,j)*katWQ;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,16,j)*katWQ;
    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Txx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,5,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(6)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,17,j)*katWT;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,17,j)*katWT;
    end
    end

    %%% Txy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,6,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(7)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,18,j)*katWT;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,18,j)*katWT;
    end
    end

    %%% Tyx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,7,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(8)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,19,j)*katWT;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,19,j)*katWT;
    end
    end

    %%% Tyy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,8,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(9)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,20,j)*katWT;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,20,j)*katWT;
    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Mxx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,9,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(11)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,21,j)*katWM;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,21,j)*katWM;
    end
    end

    %%% Mxy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,10,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(12)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,22,j)*katWM;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,22,j)*katWM;
    end
    end

    %%% Myx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,11,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(13)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,23,j)*katWM;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,23,j)*katWM;
    end
    end

    %%% Myy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,12,j);
    c3=c3+1;
    if(isnan(tk1)==1 || set.sw2(14)=='0')
        continue;
    end
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./data2(i,k,24,j)*katWM;
    c2=c2+1;wwD(c2,1)=1./data2(i,k,24,j)*katWM;
    end
    end

    %%% PTQxx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,1,j);
    tk2=data2(i,k,2,j);
    tk3=data2(i,k,3,j);
    tk4=data2(i,k,4,j);

    te1=data2(i,k,13,j);
    te2=data2(i,k,14,j);
    te3=data2(i,k,15,j);
    te4=data2(i,k,16,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,25,j)=0;
    data2(i,k,33,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(16)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,1)*katWY;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,1)*katWY;

    data2(i,k,25,j)=PT(1,1);
    data2(i,k,33,j)=sigPT(1,1);
    end
    end

    %%% PTQxy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,1,j);
    tk2=data2(i,k,2,j);
    tk3=data2(i,k,3,j);
    tk4=data2(i,k,4,j);

    te1=data2(i,k,13,j);
    te2=data2(i,k,14,j);
    te3=data2(i,k,15,j);
    te4=data2(i,k,16,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,26,j)=0;
    data2(i,k,34,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(17)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,2)*katWY;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,2)*katWY;

    data2(i,k,26,j)=PT(1,2);
    data2(i,k,34,j)=sigPT(1,2);
    end
    end

    %%% PTQyx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,1,j);
    tk2=data2(i,k,2,j);
    tk3=data2(i,k,3,j);
    tk4=data2(i,k,4,j);

    te1=data2(i,k,13,j);
    te2=data2(i,k,14,j);
    te3=data2(i,k,15,j);
    te4=data2(i,k,16,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,27,j)=0;
    data2(i,k,35,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(18)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,1)*katWY;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,1)*katWY;

    data2(i,k,27,j)=PT(2,1);
    data2(i,k,35,j)=sigPT(2,1);
    end
    end

    %%% PTQyy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,1,j);
    tk2=data2(i,k,2,j);
    tk3=data2(i,k,3,j);
    tk4=data2(i,k,4,j);

    te1=data2(i,k,13,j);
    te2=data2(i,k,14,j);
    te3=data2(i,k,15,j);
    te4=data2(i,k,16,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,28,j)=0;
    data2(i,k,36,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(19)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,2)*katWY;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,2)*katWY;

    data2(i,k,28,j)=PT(2,2);
    data2(i,k,36,j)=sigPT(2,2);
    end
    end

    %%% PTTxx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,5,j);
    tk2=data2(i,k,6,j);
    tk3=data2(i,k,7,j);
    tk4=data2(i,k,8,j);

    te1=data2(i,k,17,j);
    te2=data2(i,k,18,j);
    te3=data2(i,k,19,j);
    te4=data2(i,k,20,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,29,j)=0;
    data2(i,k,37,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(21)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));

    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,1)*katWO;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,1)*katWO;

    data2(i,k,29,j)=PT(1,1);
    data2(i,k,37,j)=sigPT(1,1);
    end
    end

    %%% PTTxy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,5,j);
    tk2=data2(i,k,6,j);
    tk3=data2(i,k,7,j);
    tk4=data2(i,k,8,j);

    te1=data2(i,k,17,j);
    te2=data2(i,k,18,j);
    te3=data2(i,k,19,j);
    te4=data2(i,k,20,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,30,j)=0;
    data2(i,k,38,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(22)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(1,2)*katWO;
    c2=c2+1;wwD(c2,1)=1./sigPT(1,2)*katWO;

    data2(i,k,30,j)=PT(1,2);
    data2(i,k,38,j)=sigPT(1,2);
    end
    end

    %%% PTTyx
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,5,j);
    tk2=data2(i,k,6,j);
    tk3=data2(i,k,7,j);
    tk4=data2(i,k,8,j);

    te1=data2(i,k,17,j);
    te2=data2(i,k,18,j);
    te3=data2(i,k,19,j);
    te4=data2(i,k,20,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,31,j)=0;
    data2(i,k,39,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(23)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,1)*katWO;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,1)*katWO;

    data2(i,k,31,j)=PT(2,1);
    data2(i,k,39,j)=sigPT(2,1);
    end
    end

    %%% PTTyy
    for j=1:nb
    for k=1:nist
    tk1=data2(i,k,5,j);
    tk2=data2(i,k,6,j);
    tk3=data2(i,k,7,j);
    tk4=data2(i,k,8,j);

    te1=data2(i,k,17,j);
    te2=data2(i,k,18,j);
    te3=data2(i,k,19,j);
    te4=data2(i,k,20,j);

    Z=[tk1 tk2; tk3 tk4];
    Ze=[te1 te2; te3 te4];

    data2(i,k,32,j)=0;
    data2(i,k,40,j)=0;
    c3=c3+1;

    if(isnan(tk1)==1 || isnan(tk2)==1 || isnan(tk3)==1 || isnan(tk4)==1 || set.sw2(24)=='0')
    continue;
    end
    [PT,sigPT,~,~] = variancePTF(Z,Ze,rand(2,1),rand(2,1));
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=c3;iv(cc)=1;
    ww(c,1)=1./sigPT(2,2)*katWO;
    c2=c2+1;wwD(c2,1)=1./sigPT(2,2)*katWO;

    data2(i,k,32,j)=PT(2,2);
    data2(i,k,40,j)=sigPT(2,2);
    end
    end

    WL{i}=sparse(ix,iy,iv,c,nist*19+nist*(4+4+4+4+4)*nb);
    Wd{i}=spdiags(ww,0,c,c);

end
% Burası modifiye edilir yeni veri eklenince

%%%soldan çarp ele
WE.WL=WL;

%%%soldan çarp ağırlıklandır
WE.Wd=Wd;
WE.WD=spdiags(wwD,0,c2,c2);

%%% sağdan çarp
for i=1:nf
    c=0;
    cc=0;
    clear ix iy iv
    ix=[];iy=[];iv=[];

    for k=1:nist

    %%Cxx - Zxy,Zxx
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    %%Cxx - Qxx, Qxy Txx Txy
    tk3=nnz(find(isnan(data(i,k,1,:))==0));
    tk4=nnz(find(isnan(data(i,k,2,:))==0));
    tk5=nnz(find(isnan(data(i,k,5,:))==0));
    tk6=nnz(find(isnan(data(i,k,6,:))==0));
    tk7=nnz(find(isnan(data(i,k,25,:))==0));
    tk8=nnz(find(isnan(data(i,k,26,:))==0));
    tk9=nnz(find(isnan(data(i,k,37,:))==0));

    if((isnan(tk1)==1 && isnan(tk2)==1 && tk3==0 && tk4==0 && tk7==0 && tk8==0 && tk9==0)  ...
        || set.sw(5)=='0' || (set.sw(2)=='0' && set.sw(1)=='0' && set.sw2(1)=='0' && set.sw2(2)=='0' && set.sw2(6)=='0' && set.sw2(7)=='0' && set.sw(18)=='0' && set.sw(19)=='0' && set.sw(25)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+1;iv(cc)=1;
    end

    %%Cxy - Zxx,Zxy
    tk1=data(i,k,1);
    tk2=data(i,k,2);
    %%Cxy - Qxx, Qxy Txx Txy
    tk3=nnz(find(isnan(data(i,k,1,:))==0));
    tk4=nnz(find(isnan(data(i,k,2,:))==0));
    tk5=nnz(find(isnan(data(i,k,5,:))==0));
    tk6=nnz(find(isnan(data(i,k,6,:))==0));
    tk7=nnz(find(isnan(data(i,k,25,:))==0));
    tk8=nnz(find(isnan(data(i,k,26,:))==0));
    tk9=nnz(find(isnan(data(i,k,37,:))==0));


    if((isnan(tk1)==1 && isnan(tk2)==1 && tk3==0 && tk4==0 && tk7==0 && tk8==0 && tk9==0) ...
       || set.sw(5)=='0' || (set.sw(2)=='0' && set.sw(1)=='0' && set.sw2(1)=='0' && set.sw2(2)=='0' && set.sw2(6)=='0' && set.sw2(7)=='0' && set.sw(18)=='0' && set.sw(19)=='0' && set.sw(25)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+2;iv(cc)=1;
    end

    %%Cyx - Zyy,Zyx
    tk1=data(i,k,3);
    tk2=data(i,k,4);
    %%Cyx - Qyy, Qyx Tyy Tyx
    tk3=nnz(find(isnan(data(i,k,3,:))==0));
    tk4=nnz(find(isnan(data(i,k,4,:))==0));
    tk5=nnz(find(isnan(data(i,k,7,:))==0));
    tk6=nnz(find(isnan(data(i,k,8,:))==0));
    tk7=nnz(find(isnan(data(i,k,27,:))==0));
    tk8=nnz(find(isnan(data(i,k,28,:))==0));
    tk9=nnz(find(isnan(data(i,k,37,:))==0));

    if((isnan(tk1)==1 && isnan(tk2)==1 && tk3==0 && tk4==0 && tk7==0 && tk8==0 && tk9==0) ...
       || set.sw(5)=='0' || (set.sw(3)=='0' && set.sw(4)=='0'  && set.sw2(3)=='0' && set.sw2(4)=='0' && set.sw2(8)=='0' && set.sw2(9)=='0' && set.sw(20)=='0' && set.sw(21)=='0' && set.sw(25)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+3;iv(cc)=1;
    end

    %%Cyy - Zyx,Zyy
    tk1=data(i,k,3);
    tk2=data(i,k,4);
    tk3=nnz(find(isnan(data(i,k,3,:))==0));
    tk4=nnz(find(isnan(data(i,k,4,:))==0));
    tk5=nnz(find(isnan(data(i,k,7,:))==0));
    tk6=nnz(find(isnan(data(i,k,8,:))==0));
    tk7=nnz(find(isnan(data(i,k,27,:))==0));
    tk8=nnz(find(isnan(data(i,k,28,:))==0));
    tk9=nnz(find(isnan(data(i,k,37,:))==0));

    if((isnan(tk1)==1 && isnan(tk2)==1 && tk3==0 && tk4==0 && tk7==0 && tk8==0 && tk9==0) ...
        || set.sw(5)=='0' || (set.sw(3)=='0' && set.sw(4)=='0' && set.sw2(3)=='0' && set.sw2(4)=='0' && set.sw2(8)=='0' && set.sw2(9)=='0' && set.sw(20)=='0' && set.sw(21)=='0' && set.sw(25)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+4;iv(cc)=1;
    end

    %%Czx - Tzy,Tzx
    tk1=data(i,k,5);
    tk2=data(i,k,6);
    if((isnan(tk1)==1 && isnan(tk2)==1)  || set.sw(9)=='0' || (set.sw(7)=='0' && set.sw(8)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+5;iv(cc)=1;
    end

    %%Czy - Tzx
    tk1=data(i,k,5);
    tk2=data(i,k,6);
    if((isnan(tk1)==1 && isnan(tk2)==1)  || set.sw(9)=='0' || (set.sw(7)=='0' && set.sw(8)=='0'))
    else
    c=c+1;
    cc=cc+1;ix(cc)=c;iy(cc)=(k-1)*6+6;iv(cc)=1;
    end
    end

    WR{i}=transpose(sparse(ix,iy,iv,c,nist*6));
end

if(set.freqdependent==0)
    ix0=[];iy0=[];
    for i=1:nf
    [ix,iy,iv]=find(transpose(WR{i}));
    ix0=[ix0;ix];
    iy0=[iy0;iy];
    end

    iy0=unique(iy0);iy0=sort(iy0);
    ix0=1:length(iy0);
    iv0=ones(size(ix0));


    for i=1:nf
        ara=max(ix0);
        if(isempty(ara)==1)
            ara=0;
        end
        WR{i}=transpose(sparse(ix0,iy0,iv0,ara,nist*6));
    end

end

D=[];
for i=1:length(base.D)
ara=base.D{1};
ara=ara';
ara=ara(:);
D=[D;ara];
end
base.D0=D;

liste=[];
if(set.freqdependent==0)
    D=base.D0;
    al=size(WR{1}',2);
    D=D(1:al);
    D=WR{1}'*D;
    base.D0=D;
else
    ek1=0;
    ek2=0;
    D=base.D0;
    DD=[];
    for i=1:nf
    ekle=size(WR{i}',2);

    liste(i,:)=zeros(1,ekle);
    koy=WR{i}'*[1:ekle]';
    liste(i,koy)=[1:length(koy)]+ek2;

    carp=D(ek1+1:ek1+ekle);
    topla=WR{i}'*carp;
    DD=[DD;topla];
    ek1=ek1+ekle;
    ek2=max(liste(i,:));
    end
    base.D0=DD;
end

base.liste=liste;

if(set.complexdist==1)
base.D0=[real(base.D0);imag(base.D0)];
end

%%%%  sağdan çarp ele distorsiyon
WE.WR=WR;

totC=0;
if(set.freqdependent==0)
    [r1]=size(WR{1},2);
    totC=totC+r1;
else
    for i=1:nf
    [r1]=size(WR{i},2);
    totC=totC+r1;
    end
end

if(set.complexdist==1)
    totC=totC*2;
end

totD=0;
d=[];
for i=1:nf

[r1]=size(WL{i},1);
totD=totD+r1;

carp=[transpose(squeeze(data(i,:,1)));...
      transpose(squeeze(data(i,:,2)));...
      transpose(squeeze(data(i,:,3)));...
      transpose(squeeze(data(i,:,4)));...
      transpose(squeeze(data(i,:,5)));...
      transpose(squeeze(data(i,:,6)));...
      transpose(squeeze(data(i,:,13)));...
      transpose(squeeze(data(i,:,14)));...
      transpose(squeeze(data(i,:,15)));...
      transpose(squeeze(data(i,:,16)));...
      transpose(squeeze(data(i,:,17)));...
      transpose(squeeze(data(i,:,18)));...
      transpose(squeeze(data(i,:,25)));...
      transpose(squeeze(data(i,:,26)));...
      transpose(squeeze(data(i,:,27)));...
      transpose(squeeze(data(i,:,28)));...
      transpose(squeeze(data(i,:,33)));...
      transpose(squeeze(data(i,:,34)));...
      transpose(squeeze(data(i,:,37)))];

for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,1,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,2,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,3,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,4,j)))];
end

for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,5,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,6,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,7,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,8,j)))];
end

for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,9,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,10,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,11,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,12,j)))];
end

for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,25,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,26,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,27,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,28,j)))];
end

for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,29,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,30,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,31,j)))];
end
for j=1:nb
carp=[carp;transpose(squeeze(data2(i,:,32,j)))];
end

ara=WL{i}*carp;
d=[d;ara];
end

base.d=d;
base.D1=base.D0;

end
