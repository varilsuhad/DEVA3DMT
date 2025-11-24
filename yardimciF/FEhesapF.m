function [Ex,Ey,Ez,Hx,Hy,Hz] = FEhesapF(EL,NK,xx,f,yuzey,set)

mu=4*pi*10^-7;


tot=size(yuzey,1);

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


for ii=1:tot
    el=yuzey(ii,1);
    
    i=yuzey(ii,2);
    j=yuzey(ii,3);
    k=yuzey(ii,4);
    a1=yuzey(ii,5);
    b1=yuzey(ii,6);   
    c1=yuzey(ii,7);


 

    [x,y,z,lx,ly,lz ] = xyzlxlylz(x,y,z,lx,ly,lz,NK,i,j,k);  

            % x(1:2)=x(1:2)/lx(1);
            % x(3:4)=x(3:4)/lx(2);
            % x(5:6)=x(5:6)/lx(3);
            % x(7:8)=x(7:8)/lx(4);
            % 
            % y([1 4])=y([1 4])/ly(1);
            % y([5 8])=y([5 8])/ly(2);
            % y([2 3])=y([2 3])/ly(3);
            % y([6 7])=y([6 7])/ly(4);    
            % 
            % z([1 5])=z([1 5])/lz(1);
            % z([2 6])=z([2 6])/lz(2);
            % z([4 8])=z([4 8])/lz(3);
            % z([3 7])=z([3 7])/lz(4);              


            % lx2=lx;
            % ly2=ly;
            % lz2=lz;   

            % %%%%%%%%%%%BÜYÜK DEĞİŞİKŞLİK
            % lx(:)=1;
            % ly(:)=1;
            % lz(:)=1;   


    [J,~] = Jabc3DMTF( x,y,z,a1,b1,c1,J);
    [R,N,M]=rotNUnF(lx,ly,lz,J,a1,b1,c1,R,N,M);
    no=EL(el,1:20);
    no=no(:);
    xx1=xx(no);
   
    acix1=atand( (z(2)-z(1))/(x(2)-x(1)));
    acix2=atand( (z(3)-z(4))/(x(3)-x(4)));
    acix=-(acix1+acix2)/2;
    
    aciy1=atand( (z(4)-z(1))/(y(4)-y(1)));
    aciy2=atand( (z(3)-z(2))/(y(3)-y(2)));
    aciy=-(aciy1+aciy2)/2;
    
    if (set.duzolcum==1)
    Rx=[1 0 0; 0 1 0; 0 0 1];    
    Ry=[1 0 0; 0 1 0; 0 0 1];    
    else
    Rx=[cosd(acix) 0 -sind(acix);...
        0 1 0;...
        sind(acix) 0 cosd(acix)];
    Ry= [1 0 0; ...
        0 cosd(aciy) -sind(aciy);...
        0 sind(aciy) cosd(aciy)];  %%x etrafında döner
    end
    
    aci=-set.olcurotasyonxy;
    Rz=[cosd(aci) -sind(aci) 0 ; ...
            sind(aci) cosd(aci) 0 ;...
        0 0 1];

    
    
    
    
    AA1=-(Rx*Ry*Rz*N')'*sqrt(-1)*2*pi;
    AA2=(Rx*Ry*Rz*M')';
    AA3=(Rz*R')'/mu;

    % AA2(1:2,1)=AA2(1:2,1)/lx2(1);
    % AA2(3:4,1)=AA2(3:4,1)/lx2(2);
    % AA2(5:6,1)=AA2(5:6,1)/lx2(3);
    % AA2(7:8,1)=AA2(7:8,1)/lx2(4);
    % 
    % AA2([1 4],2)=AA2([1 4],2)/ly2(1);
    % AA2([5 8],2)=AA2([5 8],2)/ly2(2);
    % AA2([2 3],2)=AA2([2 3],2)/ly2(3);
    % AA2([6 7],2)=AA2([6 7],2)/ly2(4); 
    % 
    % AA2([1 5],3)=AA2([1 5],3)/lz2(1);
    % AA2([2 6],3)=AA2([2 6],3)/lz2(2);
    % AA2([4 8],3)=AA2([4 8],3)/lz2(3);
    % AA2([3 7],3)=AA2([3 7],3)/lz2(4);  

    % AA2(1:2,1)=AA2(1:2,1)*lx2(1);
    % AA2(3:4,1)=AA2(3:4,1)*lx2(2);
    % AA2(5:6,1)=AA2(5:6,1)*lx2(3);
    % AA2(7:8,1)=AA2(7:8,1)*lx2(4);
    % 
    % AA2([1 4],2)=AA2([1 4],2)*ly2(1);
    % AA2([5 8],2)=AA2([5 8],2)*ly2(2);
    % AA2([2 3],2)=AA2([2 3],2)*ly2(3);
    % AA2([6 7],2)=AA2([6 7],2)*ly2(4); 
    % 
    % AA2([1 5],3)=AA2([1 5],3)*lz2(1);
    % AA2([2 6],3)=AA2([2 6],3)*lz2(2);
    % AA2([4 8],3)=AA2([4 8],3)*lz2(3);
    % AA2([3 7],3)=AA2([3 7],3)*lz2(4);              
        
    
    

    Ex(ii)=(transpose(AA1(:,1))*xx1(1:12))*f-AA2(:,1)'*xx1(13:20);
    Ey(ii)=(transpose(AA1(:,2))*xx1(1:12))*f-AA2(:,2)'*xx1(13:20);
    Ez(ii)=(transpose(AA1(:,3))*xx1(1:12))*f-AA2(:,3)'*xx1(13:20);
    
    
    Hx(ii)=AA3(:,1)'*xx1(1:12);
    Hy(ii)=AA3(:,2)'*xx1(1:12);
    Hz(ii)=AA3(:,3)'*xx1(1:12);

end


end

