% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [set] = settings3DMTF(nw,~)

set.nworkers=0;

if(nargin==1)

p = gcp('nocreate');

if (isempty(p))
    c=parpool(nw);
    else
    if(p.NumWorkers~=nw)
    delete(gcp('nocreate'));
    parpool(nw);
    end
end
set.nworkers=nw;

warning('off','all');
pctRunOnAll warning('off','all')

end

%Zxx,Zyx,Zxy,Zyy,DM - Tzx,Tzy,DMT - Sxx,Sxy,Syx,Syy,Szx,Szy -
%Pxx,Pxy,Pyx,Pyy,Pzx,Pzy - Det
set.sw='11111-000-000000-000000-0';
%Qxx,Qxy,Qyx,Qyy - Txx,Txy,Tyx,Tyy - Mxx,Mxy,Myx,Myy - Sqxx,Sqxy,Sqyx,Sqyy
%- Stxx,Stxy,Styx,Styy
set.sw2='0000-0000-0000-0000-0000';

set.baseME=0;
set.singleforwardmesh=1;

set.MergeOnCPU=0;
set.ILUT=0;
set.MatrixStacking=0;
set.triplemesh=0;
set.SPpreconditioning=0;
set.zeroinitialx=0;
set.rationalvar=0;
set.freqstacking=1;
set.GNlimit=10^-2;
set.useGPUs=1;
set.MergeOnCPU=0;

set.stagdetect=400;
set.stackmatrix=1;
set.relresForward={};
set.equalgrad=1;
set.extrapadding=0;

%%%%Çeşitli debug vs
set.skipDiscretization=0;
set.skipMapMesh=0;

set.skipCovariance=0;
set.ekblokgoster=0;
set.GaussQuadDeg=3;
set.randDD=0;   %random deltad
set.changexy=0;  %Z lerin yerini değiştir ÖNEMLİ!
set.adilpaylastirkat=1.0;
set.saveEH=0;
set.gatherveri=0;
set.forwardElek=1;    % Düz çözümde eleme yap-1  yapma-0
set.eliminatedata=0;
set.AliciIlkBlokKoy=1;
set.ortamcizmaxderin=inf;
set.ortamcizreversecb=0;
set.skipforwardOrtamkur=0;
set.sigmaort=0.01;  %%CSEM için
set.primaryEkapa=1;
set.parametreNaN=0;
set.montecarlonumber=50000;
set.initial1Dinterp=0;
set.equalgrad=1;

set.denizblok=0;  %%%% sıfır yaptım 2 idi
set.denizkarablok=0;

%%% Forward calculation
set.rotatedH=0;
set.AliciOrtada=1;  %% ikiside 0 ise alıcı neredeyse
set.AliciYuzeyde=0;
set.maxit=200;
set.limForward=10^-9;
set.limDiger=10^-9;
set.limDigerGN=10^-8;
set.gpu=1;
set.duzolcum=1;    %% Rotasyon yok elektrik alanda düşey olarak
set.olcurotasyonxy=0;  %% Elektrik ve manyetik alanı x-y düzleminde rotasyona uğrat
%%%Kapadım ki yanlışlıkla sıfır olmasın

%%% Hybrid
set.hybridpadblok=2;
set.cizdirFEFD=1;
set.allFE=0;
set.hybridfill=1;

%%% Ters çözüm
set.freqdependent=0;
set.complexdist=0;
set.preDiagAdd=0.12;
set.incompleteNE='ichol';  % 'ilu'  yada 'ichol'
set.lambda=300;
set.kappa=set.lambda*30;
set.kappaPlambda=1;

%%% Durdurma ve İlerleme
set.limitkappa=-2;
set.kappabol=2;
set.lambdabol=2;
set.stopndp=10;
set.stopkappa=15;
set.stopposrdmd=3;
set.maxitInv=200;
set.limitndp=5;   %0.1
set.stoprms=0.001;
set.stopmisfit=0.001;

%%% Başlangıç modeli
set.havaro=10^8;
set.waterro=0.33;
set.sdro=100; %En küçük blok için
set.iro=100;  %Matrisi bunla doldur
set.bro=100;  %Sınır hesabı için
set.dzks=1.1; %dz arttırma oranı
set.bounAirC=2;     %%% hava arttırma katsayısı
set.bounLandC=2;
set.bounMaxL=100;  %%% Kaç kata uzasın
set.minblok=24;  %Bölme
set.manualsd=50000;  %metre
set.autosd=1;
set.maxsd=7000000;  %metre
set.initial1D=1;
% set.skindepthortam=1;  %skin depth e göre ortam
set.meshbol=1;  %%% mesh refinement
set.meshmerge=1;
set.meshmergeincludez=0;
set.refineinterp='linear';   %%% cubic linear

%%% 1D inversion
set.MT1DInvEmp=1;  %%1D inversion 1-empedans ya da 0-log(ro)+faz
set.MT1Dfigures=1;
set.MT1Drmslim=0.2;
set.MT1Dmflim=0.5;
set.MT1Df=1000;
set.MT1Dinitalm=100;
set.MT1Dlambdalim=1;
set.MT1Drmschangestop=0.1;
set.MT1Derrorfloor=0;
set.MT1Dmaxit=15;
set.MT1Dlambda=100;

%%% Variance error floor
set.errorfloorZxy=0.03;  %%yuzde
set.errorfloorZyx=0.03;
set.errorfloorT=0.03;

%%% Süre ve yazdırmaca
set.time=0;   %%%Bireysel polarizasyon-frekans çözümü
set.relres=0;
set.frekans=0;
set.surefhazir=0;
set.sureBhazir=1;
set.sureJp=1;
set.sureJTdd=1;
set.sureForward=1; %%%Toplam  çözümü
set.sureOneIt=1; %%%Toplam  çözümü

%%%%%%% Çizdirmece %%%%%%%%%%%%%%%%%
set.errbar=[-0.05 0.05];
set.cizgisiz=1;
set.duseyImage=1;
set.yatayImage=0;
set.custombar=1;
set.cmin=1.3;
set.cmax=2.9;
set.weightedplot=1;
set.rofazcizdir=0;
set.vericizManual=1;
set.vericizstdkat=2;
set.fefdsinirciz=1;

%%%%%%%%Distortion%%%%%%%%%%%%%%%%%%
set.distortionAdd=0;
set.twist=60;
set.shear=45;
set.anisotropy=1;
set.sitegain=0;

%%%%%%%%Gürültü%%%%%%%%%%%%%%%%%%
set.noiseAdd=0;
set.noiseZ=0.02;
set.noiseT=0.03;
set.tipperrand=0;  %Distortion için

%%% Kayıt ve sayaçlar
set.rdmd=[];
set.modelr=[];
set.distr=[];
set.stop=0;
set.kappaC=0;
set.ndpC=0;
set.RMSall=[];
set.MISFIT=[];
set.lambdaIt=[];
set.kappaIt=[];
set.rs=0;
set.forwardSaveTime=[];
set.JTddSaveTime=[];
set.JpSaveTime=[];
set.OneItSaveTime=[];

% set.RMSZ=[];
% set.RMST=[];
% set.RMSpt=[];
% set.RMSpv=[];

%%% sigma katsayı
set.katWZ=1;
set.katWW=sqrt(1);
set.katWPT=sqrt(2);
set.katWPV=sqrt(2);
set.katWQ=1;
set.katWT=1;
set.katWM=sqrt(10);
set.katWY=1;
set.katWO=1;
set.katPA=sqrt(2);
set.katPB=sqrt(2);
set.katWD=sqrt(3);

set.minzok=0;

end

