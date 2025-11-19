function [FE,FD] = findDistorted3DMTF(NK,set)
%UNTÝTLED Summary of this function goes here

[nny1,nnx1,nnz1,~]=size(NK);

ny=nny1-1;
nx=nnx1-1;
nz=nnz1-1;


FE=zeros(ny,nx,nz);
FD=zeros(ny,nx,nz);


jmins=ny;
jmaxs=1;

imins=nx;
imaxs=1;

kmins=nz;
kmaxs=1;

a1=zeros(1,4);
a2=zeros(1,4);
u1=zeros(1,4);
u2=zeros(1,4);

for i=1:nx
    for j=1:ny
        for k=1:nz
            u1(1)=NK(j,i,k,3);
            u1(2)=NK(j,i+1,k,3);
            u1(3)=NK(j+1,i+1,k,3);
            u1(4)=NK(j+1,i,k,3);
            
            u2(1)=NK(j,i,k+1,3);
            u2(2)=NK(j,i+1,k+1,3);
            u2(3)=NK(j+1,i+1,k+1,3);
            u2(4)=NK(j+1,i,k+1,3);            
            
            
            f1=abs(u1(2)-u1(1))+abs(u1(3)-u1(2))+abs(u1(4)-u1(3))+abs(u1(4)-u1(1));
            f2=abs(u2(2)-u2(1))+abs(u2(3)-u2(2))+abs(u2(4)-u2(3))+abs(u2(4)-u2(1));
            
%             a1=[u1(2:end) u1(1)];
%             f1=sum(abs((u1-a1)));
%             
%             a2=[u2(2:end) u2(1)];
%             f2=sum(abs((u2-a2)));
%             

            if (f1>10^-7 || f2>10^-7)
           
            FE(j,i,k)=1;            
            
            jmaxs=max(jmaxs,j);
            jmins=min(jmins,j);
            
            imaxs=max(imaxs,i);
            imins=min(imins,i); 
            
            kmaxs=max(kmaxs,k);
            kmins=min(kmins,k);            
            
            end
        end
    end
end


if(isfield(set,'hybridfill')==0)
    set.hybridfill=1;
end

if(set.hybridfill==1)

for k=1:nz 
    aa=squeeze(FE(:,:,k));
    for j=1:ny
       ind1=find(aa(j,:)==1);
       ind2=min(ind1):max(ind1);
       aa(j,ind2)=2;
       aa(j,ind1)=1;      
       
    end
    for i=1:nx
       ind1=find(aa(:,i)==1);
       ind2=min(ind1):max(ind1);
       aa(ind2,i)=2;
       aa(ind1,i)=1;        
    end   
    FE(:,:,k)=aa;
end


for j=1:ny
    aa=squeeze(FE(j,:,:));
    for i=1:nx
       ind1=find(aa(i,:)==1);
       ind2=min(ind1):max(ind1);
       aa(i,ind2)=2;
       aa(i,ind1)=1;       
    end
    
    for k=1:nz
       ind1=find(aa(:,k)==1);
       ind2=min(ind1):max(ind1);
       aa(ind2,k)=2;
       aa(ind1,k)=1;         
    end  
    FE(j,:,:)=aa;
end

for i=1:nx
    aa=squeeze(FE(:,i,:));
    for j=1:ny
       ind1=find(aa(j,:)==1);
       ind2=min(ind1):max(ind1);
       aa(j,ind2)=2;
       aa(j,ind1)=1;
    end
    for k=1:nz
       ind1=find(aa(:,k)==1);
       ind2=min(ind1):max(ind1);
       aa(ind2,k)=2;
       aa(ind1,k)=1;       
    end   
    FE(:,i,:)=aa;
end
end
ek=set.hybridpadblok;


% 
% FE2=FE;
% if(ek>0)
% for i=1:nx
%     for j=1:ny
%         for k=1:nz
%             
%          bak=FE(j,i,k);
%          if(bak==0)
%              continue;
%          end
%          
%          j1=j-ek;
%          j2=j+ek;
%          
%          i1=i-ek;
%          i2=i+ek;      
%          
%          k1=k-ek;
%          k2=k+ek;   
%          
%          if(j1<1);j1=1;end
%          if(i1<1);i1=1;end
%          if(k1<1);k1=1;end
%          if(j2>ny);j2=ny;end
%          if(i2>nx);i2=nx;end
%          if(k2>nz);k2=nz;end         
%          
%          al=FE(j1:j2,i1:i2,k1:k2);
%          al(al==0)=3;
%          FE2(j1:j2,i1:i2,k1:k2)=al;
% 
%         end
%     end
% end
% end
% FE=FE2;


ind=find(FE~=0);

[iy,ix,iz]=ind2sub(size(FE),ind);

for i=-ek:ek
    for j=-ek:ek
        for k=-ek:ek

        ix1=ix-i;
        iy1=iy-j;
        iz1=iz-k;

        ind=(ix1<1 | ix1>nx | iy1<1 | iy1>ny | iz1<1 | iz1>nz);

        ix1(ind)=[];
        iy1(ind)=[];
        iz1(ind)=[];

        ind=sub2ind(size(FE),iy1,ix1,iz1);
        ind1=FE(ind)==0;
        FE(ind(ind1))=3;

        end
    end
end




ind=FE==0;
FD(ind)=1;



end

