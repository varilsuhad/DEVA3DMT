% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Build the core mesh, receiver, and frequency structures needed for forward and inverse runs.
function [base,set] = basekurSF(x,y,z,f,recv,set,data,data2)

if(exist('data')==1 && exist('data2')==0)
[nf,ns,~]=size(data);
data2=NaN(nf,ns,24,0);
end

st=tic;

% [x,y,z] = refinemesh3DF(x,y,z,set);
recv(:,3)=-(recv(:,3));
dx=(x(2:end)-x(1:end-1));
dy=(y(2:end)-y(1:end-1));

if(size(recv,2)==3)
al=recv(:,1:3);
recv=zeros(size(al,1),4);
recv(:,1:3)=al;
end

%     [data,f,dz,zmax,roi,roa] = anablokDataF(data,f,set);
% set.minblok=8;  %%%EKLENDI
[data,f,dz,zmax,roi,roa] = anablokDataF(data,f,set);
ny=length(dy);
nx=length(dx);
nz=length(dz);
ro=ones(ny,nx,nz);
for i=1:nx
    for j=1:ny
    ro(j,i,:)=roi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nro,nhx,nhy,nhz,ekblokx,ekbloky,ebhava,ekblokz] = blokayar3DF(ro,dx,dy,dz,min(f),set);

[ro,dz,ekblokz] = moreParameterF(ro,ebhava,ekblokz,dz,nhz);  %Yeni Eklendi

[EL,NK]=ELNK3DMTF(nro,nhx,nhy,nhz);
xcorr=sum(nhx(1:ekblokx))-x(1);
ycorr=sum(nhy(1:ekbloky))-y(1);
zcorr=sum(nhz(1:ebhava));

% [NK,buyukz1,buyukz2]=topoEkleW3DMTF(NK,ebhava,z,nhz,0.5);    %%%%2

[NK,buyukz]=topoEkle3DMTF(NK,ebhava,z,nhz,0.25);  %%%%1
NK(:,:,:,1)=NK(:,:,:,1)-xcorr;
NK(:,:,:,2)=NK(:,:,:,2)-ycorr;
NK(:,:,:,3)=NK(:,:,:,3)-zcorr;

[NK,ro] = watersurfacesetF(NK,ebhava,buyukz,ro,set);

%%%%%%%%%%%%%%%%%%m buradan
[EL,param,m,params,Re]=parametre3DMTF(EL,NK,ro,ekblokx,ekbloky,ekblokz,ebhava); %m direk sigma

% [NK,EL] = waterparameterF(EL,NK,buyukz1,buyukz2,ebhava);  %%%%2
% pause;

[FE,FD] = findDistorted3DMTF(NK,set);

% % if (set.skipCovariance~=1)
[C]=kovaryans3DMTF(NK,param,ekblokx,ekbloky,ekblokz,ebhava);  %%% Şunu sadeleştir
% else
% C=0;
% end

% if(set.allFE==1)
% FE=ones(size(FE));
% FD=zeros(size(FD));
% end

% o1=nnz(FD);
% o2=nnz(FE);
% o3=o1+o2;
% fprintf('FE=%%%.2f and FD=%%%.2f regions\n',o2/o3*100,o1/o3*100);


[yuzey]=yuzeybul3DMTF(recv,NK,EL,ebhava,set);


base.ro=ro;

base.m=m;%%%%logunu almayı unutma kovaryansla kullanırken
base.yuzey=yuzey;
base.EL=EL;
base.NK=NK;
base.nro=nro;

base.Adx=0;
base.Ady=0;

% L=L3DMTF(base.EL,base.NK,base.yuzey,set);
[DM,D]=CcreateHYB3DF(yuzey,f);

base.C1=C;
base.DM=DM;
base.D=D;
base.totP=length(Re*m);
base.totP0=length(m);

base.f=f;

ind=find(recv(:,4)~=0);
[~,ss]=sort(recv(ind,4));
ind=ind(ss);
base.brecvlist=ind;
% base.WM=speye(length(m));

base.WE=[];

[WE,totD,totC,base]=elimination3DMTALLFFF(data,data2,set,base);

base.totC=totC;
base.totD=totD;

base.WE=WE;

fprintf('Parameters(m)=%d, Distortion(C)=%d, Complex?=%d ->in Total=%d\n',length(Re*m),totC,set.complexdist,length(Re*m)+totC);
fprintf('Data=%d, frequencies=%d, stations=%d \n',totD,size(data,1),size(data,2));

base.C2=spdiags(ones(totC,1),0,totC,totC);
zero=sparse(size(base.C1,1),totC);

base.CC=[set.lambda*base.C1'*base.C1 zero; zero' set.kappa*base.C2];
base.CS=[sqrt(set.lambda)*base.C1 zero; zero' sqrt(set.kappa)*base.C2];
ek1=speye(base.totP)*set.preDiagAdd;
ek2=speye(base.totC)*set.preDiagAdd;

base.M=[set.lambda*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)];

if(set.baseME==1)
base.ME=base.M*base.M/10^100;
else
base.ME=base.M/10^100;
end
base.M=base.ME+base.M;

U=ichol(base.M);
base.ML=U;
base.MU=transpose(U);

base.ekblokx=ekblokx;
base.ekbloky=ekbloky;
base.ebhava=ebhava;

[nyo,nxo,nzo]=size(base.ro);
[Sr,Sc,Sv,N]=meshmergeB(int32(nxo),int32(nyo),int32(nzo),int32(set.meshmerge));
WM=csr2sparse(gather(Sv),gather(Sr),gather(Sc),N);
clear Sr Sc Sv

base.recv=recv;

%%%% Diğer %%%%
base.dx=dx(:);
base.dy=dy(:);
base.dz=dz(:);
base.FE=FE;
base.FD=FD;
base.smooth=WM;
base.nhx=nhx;
base.nhy=nhy;
base.nhz=nhz;

% [pay] = ortampaylastirF(f,set);

d=base.d;
brecvlist=base.brecvlist;

o1=nnz(FD);
o2=nnz(FE);
o3=o1+o2;
fprintf('Regions -> FE=%%%.2f and FD=%%%.2f\n',o2/o3*100,o1/o3*100);

if((o2/o3*100)>80)
FE=ones(size(FE));
FD=zeros(size(FD));
fprintf('Converted to Regions -> FE=%%%.2f and FD=%%%.2f\n',100,0);
end

pause(0.1);

nx=length(nhx);
ny=length(nhy);
nz=length(nhz);
NKg=gpuArray(NK);
ELg=gpuArray(EL);
FDg=gpuArray(FD);
FEg=gpuArray(FE);

nx=length(nhx);
ny=length(nhy);
nz=length(nhz);
NKg=gpuArray(NK);
ELg=gpuArray(EL);
FDg=gpuArray(FD);
FEg=gpuArray(FE);

[base.AK1v,base.AK1c,base.AK1r,base.AK2v,base.AK2c,base.AK2r,base.AK3v,base.AK3c,base.AK3r,base.AK0v,base.AK0c,base.AK0r,...
base.AM1v,base.AM1c,base.AM1r,base.AM2v,base.AM2c,base.AM2r,base.AM3v,base.AM3c,base.AM3r,base.AM4v,base.AM4c,base.AM4r, ...
base.AD1v,base.AD1c,base.AD1r,base.AD2v,base.AD2c,base.AD2r,base.AD3v,base.AD3c,base.AD3r,base.AD4v,base.AD4c,base.AD4r,...
base.AD5v,base.AD5c,base.AD5r,base.AD6v,base.AD6c,base.AD6r,base.AD7v,base.AD7c,base.AD7r,base.AD8v,base.AD8c,base.AD8r,...
base.AL1v,base.AL1c,base.AL1r,base.AL2v,base.AL2c,base.AL2r,base.AL3v,base.AL3c,base.AL3r,base.AL4v,base.AL4c,base.AL4r,...
base.AL5v,base.AL5c,base.AL5r,base.AL6v,base.AL6c,base.AL6r,base.AL7v,base.AL7c,base.AL7r,base.AL8v,base.AL8c,base.AL8r,...
base.BK1v,base.BK1c,base.BK1r,base.BK2v,base.BK2c,base.BK2r,base.BK3v,base.BK3c,base.BK3r,base.BK0v,base.BK0c,base.BK0r,...
base.BM1v,base.BM1c,base.BM1r,base.BM2v,base.BM2c,base.BM2r,base.BM3v,base.BM3c,base.BM3r,base.BM4v,base.BM4c,base.BM4r,...
base.BD1v,base.BD1c,base.BD1r,base.BD2v,base.BD2c,base.BD2r,base.BD3v,base.BD3c,base.BD3r,base.BD4v,base.BD4c,base.BD4r,...
base.BD5v,base.BD5c,base.BD5r,base.BD6v,base.BD6c,base.BD6r, base.BD7v,base.BD7c,base.BD7r,base.BD8v,base.BD8c,base.BD8r,...
base.BL1v,base.BL1c,base.BL1r,base.BL2v,base.BL2c,base.BL2r,base.BL3v,base.BL3c,base.BL3r,base.BL4v,base.BL4c,base.BL4r,...
base.BL5v,base.BL5c,base.BL5r,base.BL6v,base.BL6c,base.BL6r,base.BL7v,base.BL7c,base.BL7r,base.BL8v,base.BL8c,base.BL8r,...
base.W1v,base.W1c,base.W1r,base.W2v,base.W2c,base.W2r,base.W3v,base.W3c,base.W3r,base.W4v,base.W4c,base.W4r,...
base.W5v,base.W5c,base.W5r,base.W6v,base.W6c,base.W6r,base.W7v,base.W7c,base.W7r,base.W8v,base.W8c,base.W8r, ...
base.b]=stiffnessmultiCF(NKg,int32(ELg),int32(FDg),int32(FEg),int32(nx),int32(ny),int32(nz));

base.EL=EL;
base.NK=NK;

L=L3DMTF(base.EL,base.NK,base.yuzey,set);

[DM,~]=CcreateHYB3DF(yuzey,f);

base.L=L;
% base.ekblokz=ekblokz;
% base.flist=flist;
base.flist=f;

[nyo,nxo,nzo]=size(ro);

ko=0:N-1;
base.Wc=gpuArray(int32(ko));
ko=0:N;
base.Wr=gpuArray(int32(ko));
ko=ones(N,1);
base.Wv=gpuArray(ko);

merge=set.meshmerge;
[base.WMr,base.WMc,base.WMv,base.Sr,base.Sc,base.Sv,N2]=meshmergeCF(int32(nxo),int32(nyo),int32(nzo),int32(merge),base.Wr,base.Wc,base.Wv);
base.WM=csr2sparse(gather(base.WMv),gather(base.WMr),gather(base.WMc),N);
base.Sr=[];base.Sc=[];base.Sv=[];

base.Re=Re;

[base.L1xrow,base.L1xcol,base.L1xval] = csrgpu1(base.L.L1x);
[base.L1yrow,base.L1ycol,base.L1yval] = csrgpu1(base.L.L1y);
[base.L2xrow,base.L2xcol,base.L2xval] = csrgpu1(base.L.L2x);
[base.L2yrow,base.L2ycol,base.L2yval] = csrgpu1(base.L.L2y);
[base.Hxrow,base.Hxcol,base.Hxval] = csrgpu1(base.L.L3x);
[base.Hyrow,base.Hycol,base.Hyval] = csrgpu1(base.L.L3y);
[base.Hzrow,base.Hzcol,base.Hzval] = csrgpu1(base.L.L3z);

[base.C1row,base.C1col,base.C1val] = csrgpu1(base.DM.C1);
[base.C2row,base.C2col,base.C2val] = csrgpu1(base.DM.C2);
[base.C3row,base.C3col,base.C3val] = csrgpu1(base.DM.C3);
[base.C4row,base.C4col,base.C4val] = csrgpu1(base.DM.C4);
[base.C5row,base.C5col,base.C5val] = csrgpu1(base.DM.C5);
[base.C6row,base.C6col,base.C6val] = csrgpu1(base.DM.C6);
base.CC=size(base.DM.C1,2);

for i=1:length(base.WE.WL)
[base.WLrow{i},base.WLcol{i},base.WLval{i}] = csrgpu1(base.WE.WL{i});
base.WLC{i}=size(base.WE.WL{i},2);
[base.WRrow{i},base.WRcol{i},base.WRval{i}] = csrgpu1(base.WE.WR{i});
base.WRC{i}=size(base.WE.WR{i},2);
end

for i=1:length(base.WE.WR)
WRsize(i)=size(base.WE.WR{i},2);
end
base.WRsize=WRsize;

N=size(base.b,1);

xfor=cell(length(base.flist),1);
xJT=cell(length(base.flist),1);
xJp=cell(length(base.flist),1);
xSoUs=cell(length(base.flist),1);
xSoUs2=cell(length(base.flist),1);
xSaUs=cell(length(base.flist),1);
ii=find(base.flist~=0);
for i=1:length(ii)
xfor{ii(i)}=zeros(2*N,1,'gpuArray');
xJT{ii(i)}=zeros(2*N,1,'gpuArray');
xJp{ii(i)}=zeros(2*N,1,'gpuArray');
xSoUs{ii(i)}=zeros(2*N,1,'gpuArray');
xSoUs2{ii(i)}=zeros(2*N,1,'gpuArray');
xSaUs{ii(i)}=zeros(2*N,1,'gpuArray');
end
base.xfor=xfor;
base.xJT=xJT;
base.xJp=xJp;
base.xSoUs=xSoUs;
base.xSoUs2=xSoUs2;
base.xSaUs=xSaUs;

base.N=size(base.L.L1x,2);

st=toc(st);
fprintf('auxiliary matrices are formed in %.2f secs\n',st);



end

