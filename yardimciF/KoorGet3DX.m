function [x1,y1,z1] = KoorGet3DX(x,y,z,a,b,c)


ka=[-1 +1 +1 -1 -1 +1 +1 -1];
kb=[-1 -1 +1 +1 -1 -1 +1 +1];
kc=[-1 -1 -1 -1 +1 +1 +1 +1];


x1=0;
y1=0;
z1=0;


for i=1:8
    
x1=x1+1/8*(1+ka(i)*a)*(1+kb(i)*b)*(1+kc(i)*c)*x(i);        
y1=y1+1/8*(1+ka(i)*a)*(1+kb(i)*b)*(1+kc(i)*c)*y(i);        
z1=z1+1/8*(1+ka(i)*a)*(1+kb(i)*b)*(1+kc(i)*c)*z(i);        
    
end



end

