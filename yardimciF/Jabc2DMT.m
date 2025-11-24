% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [JJ,dd] = Jabc2DMT( x,z,a,b,J,ok)


ka=[-1 +1 +1 -1];
kb=[-1 -1 +1 +1];

s1=0;s2=0;
t1=0;t2=0;
for i=1:4
    
    kat=ka(i)*(1+kb(i)*b);
    s1=s1+kat*x(i);s2=s2+kat*z(i);

    kat=(1+ka(i)*a)*kb(i);
    t1=t1+kat*x(i);t2=t2+kat*z(i);    
end

J(1,:)=[s1 s2 ];
J(2,:)=[t1 t2];
J=(1/4)*J;


if(ok==1)
JJ=J;
dd=det(J);    
else
JJ=inv(J);
dd=det(J);
end


end

