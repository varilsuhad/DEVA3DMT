function [ nro,nhx,nhy,nhz,ekblokx,ekbloky,ebhava,ekblokz] = blokayar3DF( ro,hx,hy,hz,minf,set)


if(isfield(set,'extrapadding')==0)
set.extrapadding=0;
end

kkkh=set.bounAirC;
kkk=set.bounLandC;
ee=set.bounMaxL;


% roort=set.bro;
roort=max(ro(:));

rohava=set.havaro;


sd=503.2*sqrt(roort/minf);


nx=length(hx);
ny=length(hy);
nz=length(hz);


lkk=log10(kkk);
lkkh=log10(kkkh);

katx=ee;
katz=ee;

%ek bloklar

xyon=sd*katx;
yyon=sd*katx;
zyon=sd*katz;

%%%
c1=0;
tot=0;
xl=log10(hx(1));
while (tot<xyon)
    c1=c1+1;
    tsx(c1)=power(10,xl+lkk*(c1-1));
    tot=sum(tsx);
end
ekblokx=c1;

c2=0;
tot=0;
yl=log10(hy(1));
while (tot<yyon)
    c2=c2+1;
    tsy(c2)=power(10,yl+lkk*(c2-1));
    tot=sum(tsy);
end
ekbloky=c2;

c3=0;
tot=0;
zl=log10(hz(end));
while (tot<zyon)
    c3=c3+1;
    tsz(c3)=power(10,zl+lkk*(c3-1));
    tot=sum(tsz);
end
ekblokz=c3;


c4=0;
tot=0;
zh=log10(hz(1));
while (tot<zyon)
    c4=c4+1;
    tsh(c4)=power(10,zh+lkkh*(c4-1));
    tot=sum(tsh);
end

ebhava=c4;


ekblokx=ekblokx+set.extrapadding;
tsx(end+1:end+set.extrapadding)=tsx(end);

ekbloky=ekbloky+set.extrapadding;
tsy(end+1:end+set.extrapadding)=tsy(end);


ekblokz=ekblokz+set.extrapadding;
tsz(end+1:end+set.extrapadding)=tsz(end);

ebhava=ebhava+set.extrapadding;
tsh(end+1:end+set.extrapadding)=tsh(end);

%%%

nhx=zeros(ekblokx*2+nx,1);
nhy=zeros(ekbloky*2+ny,1);
nhz=zeros(ekblokz+nz+ebhava,1);

for i=1:ekblokx*2+nx
    if(i<=ekblokx)
    nhx(i)=tsx(end-i+1);
    end
    if(i>ekblokx && i<=ekblokx+nx)
    nhx(i)=hx(i-ekblokx);    
    end
    
    if(i>ekblokx+nx)
    nhx(i)=tsx(i-ekblokx-nx);    
    end    
end

for i=1:ekbloky*2+ny
    if(i<=ekbloky)
    nhy(i)=tsy(end-i+1);
    end
    if(i>ekbloky && i<=ekbloky+ny)
    nhy(i)=hy(i-ekbloky);    
    end
    
    if(i>ekbloky+ny)
    nhy(i)=tsy(i-ekbloky-ny);    
    end    
end

for i=1:ekblokz+nz+ebhava
    if(i<=ebhava)
        nhz(i)=tsh(end+1-i);
    end
    if(i>ebhava && i<=ebhava+nz)
    nhz(i)=hz(i-ebhava);    
    end
    
    if(i>ebhava+nz)
    nhz(i)=tsz(i-ebhava-nz);    
    end    
end

nnx=nx+2*ekblokx;
nny=ny+2*ekbloky;
noz=nz+ekblokz+ebhava;
nro=zeros(nny,nnx,noz);

nro(:,:,1:ebhava)=rohava;  %ustu hava yap


%ortaya alan? koy
for k=1:nz
    for j=1:ny
    for i=1:nx
        
        indx=ekblokx+i;
        indy=ekbloky+j;
        nro(indy,indx,ebhava+k)=ro(j,i,k);    
        
    end 
    end 
end

%yanlara dogru
    %sag-sol
    for jj=1:ny
        for kk=1:nz
       nro(ekbloky+jj,1:ekblokx,ebhava+kk)=nro(ekbloky+jj,ekblokx+1,ebhava+kk);    
       nro(ekbloky+jj,ekblokx+nx+1:end,ebhava+kk)=nro(ekbloky+jj,ekblokx+nx,ebhava+kk);    
        end 
    end
    
    %on-arka
 
    for ii=1:nx
        for kk=1:nz
       nro(1:ekbloky,ekblokx+ii,ebhava+kk)=nro(ekbloky+1,ekblokx+ii,ebhava+kk);    
       nro(ekbloky+ny+1:end,ekblokx+ii,ebhava+kk)=nro(ekbloky+ny,ekblokx+ii,ebhava+kk);    
        end 
    end 
    
    %alt
    for ii=1:nx
        for jj=1:ny
       nro(ekbloky+jj,ekblokx+ii,ebhava+nz+1:end)=nro(ekbloky+jj,ekblokx+ii,ebhava+nz);    
        end 
    end
    
    %kat kat
    
    for k=1:nz  
    nro(1:ekbloky,1:ekblokx,ebhava+k)=ro(1,1,k);    
    nro(ekbloky+ny+1:end,1:ekblokx,ebhava+k)=ro(ny,1,k); 
    nro(1:ekbloky,ekblokx+nx+1:end,ebhava+k)=ro(1,nx,k);    
    nro(ekbloky+ny+1:end,ekblokx+nx+1:end,ebhava+k)=ro(ny,nx,k);    
    end
    
    % büyük köseler
    for k=1:ekblokz
    % sol arka   
    nro(1:ekbloky+1,1:ekblokx+1,ebhava+nz+k)=ro(1,1,nz);   
    
    %sol on
    nro(ekbloky+ny:end,1:ekblokx+1,ebhava+nz+k)=ro(ny,1,nz);   

    %sag arka
     nro(1:ekbloky+1,ekblokx+nx:end,ebhava+nz+k)=ro(1,nx,nz);   
    %sag on
     nro(ekbloky+ny:end,ekblokx+nx:end,ebhava+nz+k)=ro(ny,nx,nz);   

    end   
        
%     nro
    %alt slice
    
    % nx yönü
    for i=2:nx-1
        nro(1:ekbloky,ekblokx+i,ebhava+nz+1:end)=ro(1,i,nz);
        nro(ekbloky+ny+1:end,ekblokx+i,ebhava+nz+1:end)=ro(ny,i,nz);
    end
    
    %ny yönü
    
    for j=2:ny-1
         nro(ekbloky+j,1:ekblokx,ebhava+nz+1:end)=ro(j,1,nz);
         nro(ekbloky+j,ekblokx+nx+1:end,ebhava+nz+1:end)=ro(j,nx,nz);

    end
    

% nro

fprintf('x-dir blocks=%d\n',ekblokx);
fprintf('y-dir blocks=%d\n',ekbloky);
fprintf('Air blocks=%d base blocks=%d\n',ebhava,ekblokz);
fprintf('Total Cells=%d\n',numel(nro));
fprintf('nx=%d ny=%d nz=%d\n',size(ro,2),size(ro,1),size(ro,3));
fprintf('Padded nx=%d ny=%d nz=%d\n',size(nro,2),size(nro,1),size(nro,3));

end

