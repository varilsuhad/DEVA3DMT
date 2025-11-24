function [EL,nro,m,ro,Re] = parametre3DMTF(EL,NK,ro,ekblokx,ekbloky,ekblokz,ebhava)

[ny,nx,nz]=size(ro);


c=0;c2=0;
for i=1:nx
    for j=1:ny
        for k=1:nz
        c=c+1;

        
        if(ro(j,i,k)==-2)


        m(c,1)=3.3;
        ro(j,i,k)=-c;        
        else
        m(c,1)=1./ro(j,i,k);
        ro(j,i,k)=c;  
        c2=c2+1;
        ind(c2)=c;        
        end

        end
    end
end

Re=speye(c,c);
Re=Re(ind,:);


rohava=-1;

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
    
    fprintf('m Parameter number=%d\n',length(m));

%     c=0;
    for i=1:size(EL,1)
        ii=EL(i,21);
        jj=EL(i,22);
        kk=EL(i,23);


        if(nro(jj,ii,kk)<0 && kk>ebhava)
        EL(i,24)=abs(nro(jj,ii,kk));            
        else
        EL(i,24)=nro(jj,ii,kk);            
        end
        
%         if(EL(i,24)==-1)
%         al=NK(jj:jj+1,ii:ii+1,kk:kk+1,3);al=al(:);
%         if(mean(al)>0)
%         c=c+1;    
%         EL(i,24)=-2; 
%         nro(jj,ii,kk)=-2;        
%         end
%         end
    end


    al=nro(:,:,ebhava+1:end);
    ii=al(:)<0;
    al(ii)=-2;
    nro(:,:,ebhava+1:end)=al;
    fprintf('There are %d water blocks\n',nnz(ii));



    % ii=nro<-1;
    % nro(ii)=-2;
    % fprintf('There are %d water blocks\n',nnz(ii));
    

end

