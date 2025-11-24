% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [yuzey] = yuzeybul3DMTF(recv,NK,EL,ebhava,set)

ok=0;
x=zeros(8,1);
y=zeros(8,1);
z=zeros(8,1);
lx=zeros(4,1);
ly=zeros(4,1);
lz=zeros(4,1);
yuzey=zeros(size(recv,1),12);

for j=1:size(EL,1)
    para=EL(j,24);
    kk=EL(j,23);

    if(para==-1)
        continue;
    end

    if (kk~=ebhava+1)
        continue;
    end

    ii=EL(j,21);
    jj=EL(j,22);

    [x,y,z,~,~,~] = xyzlxlylz(x,y,z,lx,ly,lz,NK,ii,jj,kk);

    mix=min(x);
    mxx=max(x);
    miy=min(y);
    mxy=max(y);
    miz=min(z);
    mxz=max(z);

    for i=1:size(recv,1)
        xb=recv(i,1);
        yb=recv(i,2);
        zb=recv(i,3);

        if (set.AliciIlkBlokKoy==1)
        zb=mean(z);
        end

        if (~(zb>=miz && zb<=mxz && xb>=mix && xb<=mxx && yb>=miy && yb<=mxy ))
        continue;
        end

        if (set.AliciIlkBlokKoy==1)
        [a,b,err] = KoorGet2DMT(x,y,xb,yb);
        c=0;
        else
        [a,b,c,err] = KoorGet3DMTF(x,y,z,xb,yb,zb);
        end

        if(err~=0)
        continue;
        end

        if (set.AliciOrtada==1)
        c=0;
        else
        c=-1;

        end


        yuzey(i,1)=j;  %eleman no
        yuzey(i,2)=ii;
        yuzey(i,3)=jj;
        yuzey(i,4)=kk;

        yuzey(i,5)=0;
        yuzey(i,6)=0;

        yuzey(i,7)=c;
        yuzey(i,8)=para;
        yuzey(i,9)=recv(i,1);
        yuzey(i,10)=recv(i,2);
        yuzey(i,11)=zb;

        if(kk~=ebhava+1)
        fprintf('Bir nokta için yerin altında yer bulundu (İlk blok değil) j\n');
        end
    end
end

for i=1:size(recv,1)
    if (yuzey(i,1)==0)
    fprintf('Bir nokta için yer bulunamadı (Hava aranmadı) i=%d \n',i);
    end
end



end

