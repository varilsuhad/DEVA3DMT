% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Assemble element connectivity and nodal coordinates for the hexahedral mesh description.
function [EL,NK] = ELNK3DMTF(ro,dx,dy,dz )

[ny,nx,nz]=size(ro);

tot=ny*nx*nz;

EL=zeros(tot,21);

for i=1:nx+1
    xx=sum(dx(1:i-1));
    for j=1:ny+1
        yy=sum(dy(1:j-1));

        for k=1:nz+1
            zz=sum(dz(1:k-1));

            NK(j,i,k,1)=xx;
            NK(j,i,k,2)=yy;
            NK(j,i,k,3)=zz;

        end
    end
end

c=0;
for i=1:nx
    for j=1:ny
        for k=1:nz

            %Ex
            if (j==1 || k==1)
            ne(1)=-1;
            else
            ne(1)=neredeAVC(i,j,k,nx,ny,nz,1);
            end

            if (k==1 || j+1==ny+1)
            ne(2)=-1;
            else
            ne(2)=neredeAVC(i,j+1,k,nx,ny,nz,1);
            end

            if(j==1 || k+1==nz+1)
            ne(3)=-1;
            else
            ne(3)=neredeAVC(i,j,k+1,nx,ny,nz,1);
            end

            if(j+1==ny+1 || k+1==nz+1)
            ne(4)=-1;
            else
            ne(4)=neredeAVC(i,j+1,k+1,nx,ny,nz,1);
            end

            %Ey
            if (i==1 || k==1)
            ne(5)=-1;
            else
            ne(5)=neredeAVC(i,j,k,nx,ny,nz,2);
            end

            if(i==1 || k+1==nz+1)
            ne(6)=-1;
            else
            ne(6)=neredeAVC(i,j,k+1,nx,ny,nz,2);
            end

            if(i+1==nx+1 || k==1)
            ne(7)=-1;
            else
            ne(7)=neredeAVC(i+1,j,k,nx,ny,nz,2);
            end

            if(i+1==nx+1 || k+1==nz+1)
            ne(8)=-1;
            else
            ne(8)=neredeAVC(i+1,j,k+1,nx,ny,nz,2);
            end

            %Ez
            if(i==1 || j==1)
            ne(9)=-1;
            else
            ne(9)=neredeAVC(i,j,k,nx,ny,nz,3);
            end

            if (i+1==nx+1 || j==1)
            ne(10)=-1;
            else
            ne(10)=neredeAVC(i+1,j,k,nx,ny,nz,3);
            end

            if(i==1 || j+1==ny+1)
            ne(11)=-1;
            else
            ne(11)=neredeAVC(i,j+1,k,nx,ny,nz,3);
            end

            if(i+1==nx+1 || j+1==ny+1)
            ne(12)=-1;
            else
            ne(12)=neredeAVC(i+1,j+1,k,nx,ny,nz,3);
            end

            if (i==1 || j==1 || k==1)
            ne(13)=-1;
            else
            ne(13)=neredeAVC(i,j,k,nx,ny,nz,4);
            end

            if (i+1==nx+1 || j==1 || k==1)
            ne(14)=-1;
            else
            ne(14)=neredeAVC(i+1,j,k,nx,ny,nz,4);
            end

            if (i+1==nx+1 || j+1==ny+1 || k==1)
            ne(15)=-1;
            else
            ne(15)=neredeAVC(i+1,j+1,k,nx,ny,nz,4);
            end

            if (i==1 || j+1==ny+1 || k==1)
            ne(16)=-1;
            else
            ne(16)=neredeAVC(i,j+1,k,nx,ny,nz,4);
            end

            if (i==1 || j==1 || k+1==nz+1)
            ne(17)=-1;
            else
            ne(17)=neredeAVC(i,j,k+1,nx,ny,nz,4);
            end

            if (i+1==nx+1 || j==1 || k+1==nz+1)
            ne(18)=-1;
            else
            ne(18)=neredeAVC(i+1,j,k+1,nx,ny,nz,4);
            end

            if (i+1==nx+1 || j+1==ny+1 || k+1==nz+1)
            ne(19)=-1;
            else
            ne(19)=neredeAVC(i+1,j+1,k+1,nx,ny,nz,4);
            end

            if (i==1 || j+1==ny+1 || k+1==nz+1)
            ne(20)=-1;
            else
            ne(20)=neredeAVC(i,j+1,k+1,nx,ny,nz,4);
            end

            c=c+1;
            EL(c,1:20)=ne;
            EL(c,21:23)=[i j k];
            EL(c,24)=0;

        end
    end
end

fprintf('DoF=%d\n',max(max(EL(:,1:20))));

end

