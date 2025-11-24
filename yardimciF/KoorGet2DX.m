% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [x1,z1] = KoorGet2DX(x,z,a,c)


ka=[-1 +1 +1 -1];
kc=[-1 -1 +1 +1];


x1=0;
z1=0;


for i=1:4
    
x1=x1+1/4*(1+ka(i)*a)*(1+kc(i)*c)*x(i);        
z1=z1+1/4*(1+ka(i)*a)*(1+kc(i)*c)*z(i);        
    
end



end

