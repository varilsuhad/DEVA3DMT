function [C] = kovaryans3DMTF(NK,param,ekblokx,ekbloky,ekblokz,ebhava)

[nny,nnx,nnz,~]=size(NK);
ny=nny-1;
nx=nnx-1;
nz=nnz-1;

x=zeros(8,1);
y=zeros(8,1);
z=zeros(8,1);
lx=zeros(4,1);
ly=zeros(4,1);
lz=zeros(4,1);

% totP=max(max(max(param)));
totP=max(param(:));

cc=0;


for j=ekbloky+1:ny-ekbloky
    for i=ekblokx+1:nx-ekblokx
         for k=ebhava+1:nz-ekblokz

            [x,y,z,lx,ly,lz] = xyzlxlylz(x,y,z,lx,ly,lz,NK,i,j,k);        
           
            c=param(j,i,k);
            
            if(c==-1)
                error('kovaryansta -1 parametre\n');
            end

            if(c==-2)
                continue;
            end            
                
            %sol
            x1=x(1);x2=x(4);x3=x(5);
            y1=y(1);y2=y(4);y3=y(5);
            z1=z(1);z2=z(4);z3=z(5);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(8);x2=x(4);x3=x(5);
            y1=y(8);y2=y(4);y3=y(5);
            z1=z(8);z2=z(4);z3=z(5);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Ssol=u1+u2;
            
            %sag
            x1=x(2);x2=x(3);x3=x(6);
            y1=y(2);y2=y(3);y3=y(6);
            z1=z(2);z2=z(3);z3=z(6);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(7);x2=x(3);x3=x(6);
            y1=y(7);y2=y(3);y3=y(6);
            z1=z(7);z2=z(3);z3=z(6);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Ssag=u1+u2;      
            
            %ust
            x1=x(1);x2=x(4);x3=x(2);
            y1=y(1);y2=y(4);y3=y(2);
            z1=z(1);z2=z(4);z3=z(2);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(3);x2=x(4);x3=x(2);
            y1=y(3);y2=y(4);y3=y(2);
            z1=z(3);z2=z(4);z3=z(2);           
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Sust=u1+u2; 
            
            %alt
            x1=x(5);x2=x(6);x3=x(8);
            y1=y(5);y2=y(6);y3=y(8);
            z1=z(5);z2=z(6);z3=z(8);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(7);x2=x(6);x3=x(8);
            y1=y(7);y2=y(6);y3=y(8);
            z1=z(7);z2=z(6);z3=z(8);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Salt=u1+u2;
            
            
            %ön
            x1=x(1);x2=x(2);x3=x(5);
            y1=y(1);y2=y(2);y3=y(5);
            z1=z(1);z2=z(2);z3=z(5);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(6);x2=x(2);x3=x(5);
            y1=y(6);y2=y(2);y3=y(5);
            z1=z(6);z2=z(2);z3=z(5);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Son=u1+u2;   
            
            
            %arka
            x1=x(4);x2=x(3);x3=x(8);
            y1=y(4);y2=y(3);y3=y(8);
            z1=z(4);z2=z(3);z3=z(8);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            u1=1/2*sqrt(s1^2+s2^2+s3^2);
            
            x1=x(7);x2=x(3);x3=x(8);
            y1=y(7);y2=y(3);y3=y(8);
            z1=z(7);z2=z(3);z3=z(8);            
            s1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
            s2=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
            s3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);            
       
            u2=1/2*sqrt(s1^2+s2^2+s3^2);
            Sarka=u1+u2;               
            
            top=0;
            %i-1
            if(i==ekblokx+1 || param(j,i-1,k)<0) 
            else
            top=top-Ssol;
            end
            %i+1
            if(i==nx-ekblokx || param(j,i+1,k)<0) 
            else
            top=top-Ssag;  
            end   
            %j-1
            if(j==ekbloky+1 || param(j-1,i,k)<0) 
            else
            top=top-Son;
            end            
            %j+1
            if(j==ny-ekbloky || param(j+1,i,k)<0) 
            else
            top=top-Sarka;
            end              
            %k-1
            if(k==ebhava+1 || param(j,i,k-1)<0) 
            else
            top=top-Sust;
            end 
            %k+1
            if(k==nz-ekblokz || param(j,i,k+1)<0) 
            else
            top=top-Salt;
            end  
            
            %%%orta
            cc=cc+1;ix(cc)=c;iy(cc)=c;iv(cc)=-1;                                                                         
%             iv(cc)=top;  

            %i-1
            if(i==ekblokx+1 || param(j,i-1,k)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j,i-1,k);iv(cc)=-Ssol/top;             
            end
            
            %i+1
            if(i==nx-ekblokx || param(j,i+1,k)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j,i+1,k);iv(cc)=-Ssag/top;                         
            end   
            
            %j-1
            if(j==ekbloky+1 || param(j-1,i,k)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j-1,i,k);iv(cc)=-Son/top;                                     
            end            
    
            %j+1
            if(j==ny-ekbloky || param(j+1,i,k)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j+1,i,k);iv(cc)=-Sarka/top;                                                 
            end              
            
            %k-1
            if(k==ebhava+1 || param(j,i,k-1)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j,i,k-1);iv(cc)=-Sust/top;                                                             
            end 
            
            %k+1
            if(k==nz-ekblokz || param(j,i,k+1)<0) 
            else
            cc=cc+1;ix(cc)=c;iy(cc)=param(j,i,k+1);iv(cc)=-Salt/top;                                                                         
            end    
            
            
        end
    end
end

C=sparse(ix,iy,iv,totP,totP);  


end

