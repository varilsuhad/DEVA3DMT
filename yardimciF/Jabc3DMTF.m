% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Form the inverse Jacobian and determinant for trilinear hexahedral mapping at a given local point.
function [JJ,dd] = Jabc3DMTF( x,y,z,a,b,c,JJ)

q1=-(1-b)*(1-c);
q2=(1-b)*(1-c);
q3=(1+b)*(1-c);
q4=-(1+b)*(1-c);
q5=-(1-b)*(1+c);
q6=(1-b)*(1+c);
q7=(1+b)*(1+c);
q8=-(1+b)*(1+c);

w1=-(1-a)*(1-c);
w2=-(1+a)*(1-c);
w3=(1+a)*(1-c);
w4=(1-a)*(1-c);
w5=-(1-a)*(1+c);
w6=-(1+a)*(1+c);
w7=(1+a)*(1+c);
w8=(1-a)*(1+c);

e1=-(1-a)*(1-b);
e2=-(1+a)*(1-b);
e3=-(1+a)*(1+b);
e4=-(1-a)*(1+b);
e5=(1-a)*(1-b);
e6=(1+a)*(1-b);
e7=(1+a)*(1+b);
e8=(1-a)*(1+b);

s1=q1*x(1)+q2*x(2)+q3*x(3)+q4*x(4)+q5*x(5)+q6*x(6)+q7*x(7)+q8*x(8);
s2=q1*y(1)+q2*y(2)+q3*y(3)+q4*y(4)+q5*y(5)+q6*y(6)+q7*y(7)+q8*y(8);
s3=q1*z(1)+q2*z(2)+q3*z(3)+q4*z(4)+q5*z(5)+q6*z(6)+q7*z(7)+q8*z(8);

t1=w1*x(1)+w2*x(2)+w3*x(3)+w4*x(4)+w5*x(5)+w6*x(6)+w7*x(7)+w8*x(8);
t2=w1*y(1)+w2*y(2)+w3*y(3)+w4*y(4)+w5*y(5)+w6*y(6)+w7*y(7)+w8*y(8);
t3=w1*z(1)+w2*z(2)+w3*z(3)+w4*z(4)+w5*z(5)+w6*z(6)+w7*z(7)+w8*z(8);

r1=e1*x(1)+e2*x(2)+e3*x(3)+e4*x(4)+e5*x(5)+e6*x(6)+e7*x(7)+e8*x(8);
r2=e1*y(1)+e2*y(2)+e3*y(3)+e4*y(4)+e5*y(5)+e6*y(6)+e7*y(7)+e8*y(8);
r3=e1*z(1)+e2*z(2)+e3*z(3)+e4*z(4)+e5*z(5)+e6*z(6)+e7*z(7)+e8*z(8);

r1=0.125*r1;
r2=0.125*r2;
r3=0.125*r3;
t1=0.125*t1;
t2=0.125*t2;
t3=0.125*t3;
s1=0.125*s1;
s2=0.125*s2;
s3=0.125*s3;

dd=(s1*(t2*r3-t3*r2)-s2*(t1*r3-r1*t3)+s3*(t1*r2-t2*r1));
JJ(1,1)=(t2*r3-r2*t3)/dd;
JJ(1,2)=-(s2*r3-r2*s3)/dd;
JJ(1,3)=(s2*t3-s3*t2)/dd;
JJ(2,1)=-(t1*r3-r1*t3)/dd;
JJ(2,2)=(s1*r3-r1*s3)/dd;
JJ(2,3)=-(s1*t3-t1*s3)/dd;
JJ(3,1)=(t1*r2-r1*t2)/dd;
JJ(3,2)=-(s1*r2-s2*r1)/dd;
JJ(3,3)=(s1*t2-t1*s2)/dd;

end

