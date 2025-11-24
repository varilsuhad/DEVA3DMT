% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [x,y,z,lx,ly,lz ] = xyzlxlylz(x,y,z,lx,ly,lz,NK,i,j,k)
x(1)=NK(j,i,k,1);
x(2)=NK(j,i+1,k,1);
x(3)=NK(j+1,i+1,k,1);
x(4)=NK(j+1,i,k,1);
x(5)=NK(j,i,k+1,1);
x(6)=NK(j,i+1,k+1,1);
x(7)=NK(j+1,i+1,k+1,1);
x(8)=NK(j+1,i,k+1,1);

y(1)=NK(j,i,k,2);
y(2)=NK(j,i+1,k,2);
y(3)=NK(j+1,i+1,k,2);
y(4)=NK(j+1,i,k,2);
y(5)=NK(j,i,k+1,2);
y(6)=NK(j,i+1,k+1,2);
y(7)=NK(j+1,i+1,k+1,2);
y(8)=NK(j+1,i,k+1,2);

z(1)=NK(j,i,k,3);
z(2)=NK(j,i+1,k,3);
z(3)=NK(j+1,i+1,k,3);
z(4)=NK(j+1,i,k,3);
z(5)=NK(j,i,k+1,3);
z(6)=NK(j,i+1,k+1,3);
z(7)=NK(j+1,i+1,k+1,3);
z(8)=NK(j+1,i,k+1,3);

lx(1)=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
lx(2)=sqrt((x(3)-x(4))^2+(y(3)-y(4))^2+(z(3)-z(4))^2);
lx(3)=sqrt((x(6)-x(5))^2+(y(6)-y(5))^2+(z(6)-z(5))^2);
lx(4)=sqrt((x(7)-x(8))^2+(y(7)-y(8))^2+(z(7)-z(8))^2);

ly(1)=sqrt((x(4)-x(1))^2+(y(4)-y(1))^2+(z(4)-z(1))^2);
ly(2)=sqrt((x(8)-x(5))^2+(y(8)-y(5))^2+(z(8)-z(5))^2);
ly(3)=sqrt((x(3)-x(2))^2+(y(3)-y(2))^2+(z(3)-z(2))^2);
ly(4)=sqrt((x(7)-x(6))^2+(y(7)-y(6))^2+(z(7)-z(6))^2);

lz(1)=sqrt((x(5)-x(1))^2+(y(5)-y(1))^2+(z(5)-z(1))^2);
lz(2)=sqrt((x(6)-x(2))^2+(y(6)-y(2))^2+(z(6)-z(2))^2);
lz(3)=sqrt((x(8)-x(4))^2+(y(8)-y(4))^2+(z(8)-z(4))^2);
lz(4)=sqrt((x(7)-x(3))^2+(y(7)-y(3))^2+(z(7)-z(3))^2);

end

