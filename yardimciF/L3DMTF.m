function [L] = L3DMTF(EL,NK,yuzey,set)


fprintf('Rotation angle(xy)=%.2f (L) \n',set.olcurotasyonxy);

mu=4*pi*10^-7;
totE=max(max(EL(:,1:20)));
totL=size(yuzey,1);

x=zeros(8,1);
y=zeros(8,1);
z=zeros(8,1);
lx=zeros(4,1);
ly=zeros(4,1);
lz=zeros(4,1);
R=zeros(12,3);
N=zeros(12,3);
M=zeros(8,3);
J=zeros(3,3);

son=0;
son2=0;
for ii=1:totL
    el=yuzey(ii,1);
    
    i=yuzey(ii,2);
    j=yuzey(ii,3);
    k=yuzey(ii,4);
    a1=yuzey(ii,5);
    b1=yuzey(ii,6);   
    c1=yuzey(ii,7);

    [x,y,z,lx,ly,lz ] = xyzlxlylz(x,y,z,lx,ly,lz,NK,i,j,k);        
    [J,~] = Jabc3DMTF( x,y,z,a1,b1,c1,J);
    [R,N,M]=rotNUnF(lx,ly,lz,J,a1,b1,c1,R,N,M);
    
    no=EL(el,1:20);
      
    acix1=atand( (z(2)-z(1))/(x(2)-x(1)));
    acix2=atand( (z(3)-z(4))/(x(3)-x(4)));
    acix=-(acix1+acix2)/2;
    
    aciy1=atand( (z(4)-z(1))/(y(4)-y(1)));
    aciy2=atand( (z(3)-z(2))/(y(3)-y(2)));
    aciy=-(aciy1+aciy2)/2;
    

%     Rx=[cosd(acix) 0 -sind(acix);...
%         0 1 0;...
%         sind(acix) 0 cosd(acix)];
%     Ry= [1 0 0; ...
%         0 cosd(aciy) -sind(aciy);...
%         0 sind(aciy) cosd(aciy)];  %%x etrafýnda döner
    
    if (set.duzolcum==1)
    Rx=[1 0 0; 0 1 0; 0 0 1];    
    Ry=[1 0 0; 0 1 0; 0 0 1];    
    else
    Rx=[cosd(acix) 0 -sind(acix);...
        0 1 0;...
        sind(acix) 0 cosd(acix)];
    Ry= [1 0 0; ...
        0 cosd(aciy) -sind(aciy);...
        0 sind(aciy) cosd(aciy)];  %%x etrafýnda döner
    end
    
    aci=-set.olcurotasyonxy;
    Rz=[cosd(aci) -sind(aci) 0 ; ...
        sind(aci) cosd(aci) 0 ;...
        0 0 1];    
    
    
    
    AA1=-(Rx*Ry*Rz*N')'*sqrt(-1)*2*pi;
    AA2=(Rx*Ry*Rz*M')';
    AA3=(Rz*R')'/mu;
    
    
    ek1=son+1;ek2=son+12;son=son+12;
    xx=ones(12,1)*ii;yy=no(1:12);
    ix1(ek1:ek2)=xx;iy1(ek1:ek2)=yy;
    ie1(ek1:ek2)=AA1(:,1);ie2(ek1:ek2)=AA1(:,2);ie3(ek1:ek2)=AA1(:,3);     
    ih1(ek1:ek2)=AA3(:,1);ih2(ek1:ek2)=AA3(:,2);ih3(ek1:ek2)=AA3(:,3);     
    
    
    ek1=son2+1;ek2=son2+8;son2=son2+8;
    xx=ones(8,1)*ii;yy=no(13:20);
    ix2(ek1:ek2)=xx;iy2(ek1:ek2)=yy;
    iv1(ek1:ek2)=AA2(:,1);iv2(ek1:ek2)=AA2(:,2);iv3(ek1:ek2)=AA2(:,3);      
  
end

%%% E'nin A'sýný      f ile çarpýlýr
L1x=sparse(ix1,iy1,ie1,totL,totE);
L1y=sparse(ix1,iy1,ie2,totL,totE);
L1z=sparse(ix1,iy1,ie3,totL,totE);
%%%E'nin v'sini       - ile çarpýlýr
L2x=sparse(ix2,iy2,iv1,totL,totE);
L2y=sparse(ix2,iy2,iv2,totL,totE);
L2z=sparse(ix2,iy2,iv3,totL,totE);
%%% H'ý                sade böyle
L3x=sparse(ix1,iy1,ih1,totL,totE);
L3y=sparse(ix1,iy1,ih2,totL,totE);
L3z=sparse(ix1,iy1,ih3,totL,totE);

L.L1x=L1x;
L.L1y=L1y;
L.L1z=L1z;

L.L2x=L2x;
L.L2y=L2y;
L.L2z=L2z;

L.L3x=L3x;
L.L3y=L3y;
L.L3z=L3z;


end

