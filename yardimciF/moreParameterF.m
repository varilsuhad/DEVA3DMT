function [ro2,dz2,ekblokz2] = moreParameterF(ro,ebhava,ekblokz,dz,nhz)


dz=dz(:);
nhz=nhz(:);


kalan=2;

if(ekblokz-kalan<0)
    error('moreParameter error\n');
end

fark=ekblokz-kalan;


dz2=[dz;nhz(ebhava+length(dz)+1:end-kalan)];

[ny,nx,nz]=size(ro);

ro2=zeros(ny,nx,nz+fark);

ro2(:,:,1:nz)=ro;

for i=1:fark
ro2(:,:,nz+i)=ro(:,:,end);
end

ekblokz2=kalan;


end