function [NK,nz1] = topoEkle3DMTF(NK,ebhava,z,nhz,ek)

[nny,nnx,~,~]=size(NK);


zo=z;

ii=find(z(:)<0);
fprintf("%d points are under sea level\n",length(ii));
z(ii)=0;



meanz=mean(z(:));
z=z-meanz;

[nny2,nnx2]=size(z);

ekblokny=(nny-nny2)/2;
ekbloknx=(nnx-nnx2)/2;

nz=zeros(2*ekblokny+nny2,2*ekbloknx+nnx2);

nz(ekblokny+1:ekblokny+nny2,ekbloknx+1:ekbloknx+nnx2)=z;

nz(1:ekblokny+1,1:ekbloknx+1)=nz(ekblokny+1,ekbloknx+1);
nz(1:ekblokny+1,ekbloknx+nnx2:end)=nz(ekblokny+1,ekbloknx+nnx2);
nz(ekblokny+nny2:end,1:ekbloknx+1)=nz(ekblokny+nny2,ekbloknx+1);
nz(ekblokny+nny2:end,ekbloknx+nnx2:end)=nz(ekblokny+nny2,ekbloknx+nnx2);

for i=1:nny2
nz(ekblokny+i,1:ekbloknx)=nz(ekblokny+i,ekbloknx+1);
nz(ekblokny+i,ekbloknx+nnx2+1:end)=nz(ekblokny+i,ekbloknx+nnx2);
end

for i=1:nnx2
nz(1:ekblokny,ekbloknx+i)=nz(ekblokny+1,ekbloknx+i);
nz(ekblokny+nny2+1:end,ekbloknx+i)=nz(ekblokny+nny2,ekbloknx+i);
end

for i=1:nnx
    for j=1:nny         
     NK(j,i,:,3)=topoEkleSubF(nhz,nz(j,i),ebhava,ek);
    end
end

NK(:,:,:,3)=NK(:,:,:,3)-meanz;
nz=nz+meanz;
% z=z+meanz;

[sy1,sx1]=size(nz);
nz1=zeros(sy1,sx1);


[sy,sx]=size(zo);


for i=1:sx
    for j=1:sy


       nz1(j+ekblokny,i+ekbloknx)=zo(j,i);
       if(i==1 && j==1)
       nz1(1:ekblokny+1,1:ekbloknx+1)=zo(j,i);
       continue;       
       elseif(i==1 && j==sy)
       nz1(end-ekblokny:end,1:ekbloknx+1)=zo(j,i);
       continue;       
       elseif(i==sx && j==1)
       nz1(1:ekblokny+1,end-ekbloknx:end)=zo(j,i);
       continue;       
       elseif(i==sx && j==sy)
       nz1(end-ekblokny:end,end-ekbloknx:end)=zo(j,i);
       continue;
       elseif(i==1)
       nz1(j+ekblokny,1:ekbloknx)=zo(j,i);
       continue;
       elseif(i==sx)
           nz1(j+ekblokny,end-ekbloknx:end)=zo(j,i);
       continue;
       elseif(j==1)
           nz1(1:ekblokny,i+ekbloknx)=zo(j,i);
       continue;
       elseif(j==sy)
           nz1(end-ekblokny:end,i+ekbloknx)=zo(j,i);
       continue;       
       end


    end
end






end

