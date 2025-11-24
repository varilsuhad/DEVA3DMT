% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Merge polarization-specific forward solutions into combined prediction vectors for selected data types.
function [ bak1,bak2 ] = Pmerge3DMTHYBS3MULTF3(base,veri,ok,vec1,vec2)

% [ bak1,bak2 ]
% [ bak1,bak2,bak3,bak4,bak5,bak6,bak7,bak8 ]
%%% ok=1 P*vec
%%% ok=2 P'*vec

f=veri.f;
fd=1./f;

bak1=0;
bak2=0;

if (ok==1)
e1=(veri.e1);
e2=(veri.e2);

%%%%%%%%%%%%%%%%%% FE %%%%%%%%%%%%%%%%%%%%%%%%%

wm1=base.WM*(base.m.*vec1);
wm2=base.WM*(base.m.*vec2);

v1=e1.*(base.BW1*wm1);
v2=e1.*(base.BW2*wm1);
v3=e1.*(base.BW3*wm1);
v4=e1.*(base.BW4*wm1);
v5=e1.*(base.BW5*wm1);
v6=e1.*(base.BW6*wm1);
v7=e1.*(base.BW7*wm1);
v8=e1.*(base.BW8*wm1);

y1=e2.*(base.BW1*wm2);
y2=e2.*(base.BW2*wm2);
y3=e2.*(base.BW3*wm2);
y4=e2.*(base.BW4*wm2);
y5=e2.*(base.BW5*wm2);
y6=e2.*(base.BW6*wm2);
y7=e2.*(base.BW7*wm2);
y8=e2.*(base.BW8*wm2);

X1=(sqrt(-1)*base.BM1*f*v1+base.BD1*v1+sqrt(-1)*base.BL1*fd*v1);
X2=(sqrt(-1)*base.BM2*f*v2+base.BD2*v2+sqrt(-1)*base.BL2*fd*v2);
X3=(sqrt(-1)*base.BM3*f*v3+base.BD3*v3+sqrt(-1)*base.BL3*fd*v3);
X4=(sqrt(-1)*base.BM4*f*v4+base.BD4*v4+sqrt(-1)*base.BL4*fd*v4);
X5=(base.BD5*v5+sqrt(-1)*base.BL5*fd*v5);
X6=(base.BD6*v6+sqrt(-1)*base.BL6*fd*v6);
X7=(base.BD7*v7+sqrt(-1)*base.BL7*fd*v7);
X8=(base.BD8*v8+sqrt(-1)*base.BL8*fd*v8);

bak11=-X1...
      -X2 ...
      -X3 ...
      -X4 ...
      -X5 ...
      -X6 ...
      -X7 ...
      -X8;

X1=(sqrt(-1)*base.BM1*f*y1+base.BD1*y1+sqrt(-1)*base.BL1*fd*y1);
X2=(sqrt(-1)*base.BM2*f*y2+base.BD2*y2+sqrt(-1)*base.BL2*fd*y2);
X3=(sqrt(-1)*base.BM3*f*y3+base.BD3*y3+sqrt(-1)*base.BL3*fd*y3);
X4=(sqrt(-1)*base.BM4*f*y4+base.BD4*y4+sqrt(-1)*base.BL4*fd*y4);
X5=(base.BD5*y5+sqrt(-1)*base.BL5*fd*y5);
X6=(base.BD6*y6+sqrt(-1)*base.BL6*fd*y6);
X7=(base.BD7*y7+sqrt(-1)*base.BL7*fd*y7);
X8=(base.BD8*y8+sqrt(-1)*base.BL8*fd*y8);

bak12=-X1 ...
      -X2 ...
      -X3 ...
      -X4 ...
      -X5 ...
      -X6...
      -X7 ...
      -X8;

v1=e1.*(base.AW1*wm1);
v2=e1.*(base.AW2*wm1);
v3=e1.*(base.AW3*wm1);
v4=e1.*(base.AW4*wm1);
v5=e1.*(base.AW5*wm1);
v6=e1.*(base.AW6*wm1);
v7=e1.*(base.AW7*wm1);
v8=e1.*(base.AW8*wm1);

y1=e2.*(base.AW1*wm2);
y2=e2.*(base.AW2*wm2);
y3=e2.*(base.AW3*wm2);
y4=e2.*(base.AW4*wm2);
y5=e2.*(base.AW5*wm2);
y6=e2.*(base.AW6*wm2);
y7=e2.*(base.AW7*wm2);
y8=e2.*(base.AW8*wm2);

X1=(sqrt(-1)*base.AM1*f*v1+base.AD1*v1+sqrt(-1)*base.AL1*fd*v1);
X2=(sqrt(-1)*base.AM2*f*v2+base.AD2*v2+sqrt(-1)*base.AL2*fd*v2);
X3=(sqrt(-1)*base.AM3*f*v3+base.AD3*v3+sqrt(-1)*base.AL3*fd*v3);
X4=(sqrt(-1)*base.AM4*f*v4+base.AD4*v4+sqrt(-1)*base.AL4*fd*v4);
X5=(base.AD5*v5+sqrt(-1)*base.AL5*fd*v5);
X6=(base.AD6*v6+sqrt(-1)*base.AL6*fd*v6);
X7=(base.AD7*v7+sqrt(-1)*base.AL7*fd*v7);
X8=(base.AD8*v8+sqrt(-1)*base.AL8*fd*v8);

bak21=-X1...
      -X2 ...
      -X3 ...
      -X4 ...
      -X5 ...
      -X6 ...
      -X7 ...
      -X8;
X1=(sqrt(-1)*base.AM1*f*y1+base.AD1*y1+sqrt(-1)*base.AL1*fd*y1);
X2=(sqrt(-1)*base.AM2*f*y2+base.AD2*y2+sqrt(-1)*base.AL2*fd*y2);
X3=(sqrt(-1)*base.AM3*f*y3+base.AD3*y3+sqrt(-1)*base.AL3*fd*y3);
X4=(sqrt(-1)*base.AM4*f*y4+base.AD4*y4+sqrt(-1)*base.AL4*fd*y4);
X5=(base.AD5*y5+sqrt(-1)*base.AL5*fd*y5);
X6=(base.AD6*y6+sqrt(-1)*base.AL6*fd*y6);
X7=(base.AD7*y7+sqrt(-1)*base.AL7*fd*y7);
X8=(base.AD8*y8+sqrt(-1)*base.AL8*fd*y8);

bak22=-X1 ...
      -X2 ...
      -X3 ...
      -X4 ...
      -X5 ...
      -X6...
      -X7 ...
      -X8;


bak1=bak11+bak21;
bak2=bak12+bak22;


elseif(ok==2)


e1=(veri.e1.');
e2=(veri.e2.');

%%%%%%%%%%%%%%%%%% FE %%%%%%%%%%%%%%%%%%%%%%%%%

X1=(vec1.'*sqrt(-1)*base.BM1*f+vec1.'*base.BD1+vec1.'*sqrt(-1)*base.BL1*fd);
X2=(vec1.'*sqrt(-1)*base.BM2*f+vec1.'*base.BD2+vec1.'*sqrt(-1)*base.BL2*fd);
X3=(vec1.'*sqrt(-1)*base.BM3*f+vec1.'*base.BD3+vec1.'*sqrt(-1)*base.BL3*fd);
X4=(vec1.'*sqrt(-1)*base.BM4*f+vec1.'*base.BD4+vec1.'*sqrt(-1)*base.BL4*fd);
X5=(vec1.'*base.BD5+vec1.'*sqrt(-1)*base.BL5*fd);
X6=(vec1.'*base.BD6+vec1.'*sqrt(-1)*base.BL6*fd);
X7=(vec1.'*base.BD7+vec1.'*sqrt(-1)*base.BL7*fd);
X8=(vec1.'*base.BD8+vec1.'*sqrt(-1)*base.BL8*fd);

bak11=-(X1.*e1)*(base.BW1)...
      -(X2.*e1)*(base.BW2)...
      -(X3.*e1)*(base.BW3)...
      -(X4.*e1)*(base.BW4)...
      -(X5.*e1)*(base.BW5)...
      -(X6.*e1)*(base.BW6)...
      -(X7.*e1)*(base.BW7)...
      -(X8.*e1)*(base.BW8);

bak11=base.m.*(bak11*(base.WM)).';

X1=(vec2.'*sqrt(-1)*base.BM1*f+vec2.'*base.BD1+vec2.'*sqrt(-1)*base.BL1*fd);
X2=(vec2.'*sqrt(-1)*base.BM2*f+vec2.'*base.BD2+vec2.'*sqrt(-1)*base.BL2*fd);
X3=(vec2.'*sqrt(-1)*base.BM3*f+vec2.'*base.BD3+vec2.'*sqrt(-1)*base.BL3*fd);
X4=(vec2.'*sqrt(-1)*base.BM4*f+vec2.'*base.BD4+vec2.'*sqrt(-1)*base.BL4*fd);
X5=(vec2.'*base.BD5+vec2.'*sqrt(-1)*base.BL5*fd);
X6=(vec2.'*base.BD6+vec2.'*sqrt(-1)*base.BL6*fd);
X7=(vec2.'*base.BD7+vec2.'*sqrt(-1)*base.BL7*fd);
X8=(vec2.'*base.BD8+vec2.'*sqrt(-1)*base.BL8*fd);

bak12=-(X1.*e2)*(base.BW1)...
      -(X2.*e2)*(base.BW2)...
      -(X3.*e2)*(base.BW3)...
      -(X4.*e2)*(base.BW4)...
      -(X5.*e2)*(base.BW5)...
      -(X6.*e2)*(base.BW6)...
      -(X7.*e2)*(base.BW7)...
      -(X8.*e2)*(base.BW8);
bak12=base.m.*(bak12*(base.WM)).';

%%%%%%%%%%%%%%%%%% FD %%%%%%%%%%%%%%%%%%%%%%%%%

X1=(vec1.'*sqrt(-1)*base.AM1*f+vec1.'*base.AD1+vec1.'*sqrt(-1)*base.AL1*fd);
X2=(vec1.'*sqrt(-1)*base.AM2*f+vec1.'*base.AD2+vec1.'*sqrt(-1)*base.AL2*fd);
X3=(vec1.'*sqrt(-1)*base.AM3*f+vec1.'*base.AD3+vec1.'*sqrt(-1)*base.AL3*fd);
X4=(vec1.'*sqrt(-1)*base.AM4*f+vec1.'*base.AD4+vec1.'*sqrt(-1)*base.AL4*fd);
X5=(vec1.'*base.AD5+vec1.'*sqrt(-1)*base.AL5*fd);
X6=(vec1.'*base.AD6+vec1.'*sqrt(-1)*base.AL6*fd);
X7=(vec1.'*base.AD7+vec1.'*sqrt(-1)*base.AL7*fd);
X8=(vec1.'*base.AD8+vec1.'*sqrt(-1)*base.AL8*fd);

bak21=-(X1.*e1)*(base.AW1)...
      -(X2.*e1)*(base.AW2)...
      -(X3.*e1)*(base.AW3)...
      -(X4.*e1)*(base.AW4)...
      -(X5.*e1)*(base.AW5)...
      -(X6.*e1)*(base.AW6)...
      -(X7.*e1)*(base.AW7)...
      -(X8.*e1)*(base.AW8);
bak21=base.m.*(bak21*(base.WM)).';

X1=(vec2.'*sqrt(-1)*base.AM1*f+vec2.'*base.AD1+vec2.'*sqrt(-1)*base.AL1*fd);
X2=(vec2.'*sqrt(-1)*base.AM2*f+vec2.'*base.AD2+vec2.'*sqrt(-1)*base.AL2*fd);
X3=(vec2.'*sqrt(-1)*base.AM3*f+vec2.'*base.AD3+vec2.'*sqrt(-1)*base.AL3*fd);
X4=(vec2.'*sqrt(-1)*base.AM4*f+vec2.'*base.AD4+vec2.'*sqrt(-1)*base.AL4*fd);
X5=(vec2.'*base.AD5+vec2.'*sqrt(-1)*base.AL5*fd);
X6=(vec2.'*base.AD6+vec2.'*sqrt(-1)*base.AL6*fd);
X7=(vec2.'*base.AD7+vec2.'*sqrt(-1)*base.AL7*fd);
X8=(vec2.'*base.AD8+vec2.'*sqrt(-1)*base.AL8*fd);

bak22=-(X1.*e2)*(base.AW1)...
      -(X2.*e2)*(base.AW2)...
      -(X3.*e2)*(base.AW3)...
      -(X4.*e2)*(base.AW4)...
      -(X5.*e2)*(base.AW5)...
      -(X6.*e2)*(base.AW6)...
      -(X7.*e2)*(base.AW7)...
      -(X8.*e2)*(base.AW8);
bak22=base.m.*(bak22*(base.WM)).';


bak1=bak11+bak21;
bak2=bak12+bak22;

elseif(ok==3)
e1=(veri.e1);
e2=(veri.e2);

e1=spdiags(e1,0,length(e1),length(e1));
e2=spdiags(e2,0,length(e2),length(e2));
%%%%%%%%%%%%%%%%%% FE %%%%%%%%%%%%%%%%%%%%%%%%%

X1=(sqrt(-1)*base.BM1*f+base.BD1+sqrt(-1)*base.BL1/f);
X2=(sqrt(-1)*base.BM2*f+base.BD2+sqrt(-1)*base.BL2/f);
X3=(sqrt(-1)*base.BM3*f+base.BD3+sqrt(-1)*base.BL3/f);
X4=(sqrt(-1)*base.BM4*f+base.BD4+sqrt(-1)*base.BL4/f);
X5=(base.BD5+sqrt(-1)*base.BL5/f);
X6=(base.BD6+sqrt(-1)*base.BL6/f);
X7=(base.BD7+sqrt(-1)*base.BL7/f);
X8=(base.BD8+sqrt(-1)*base.BL8/f);

bak11=-X1*(e1*(base.BW1*(base.WM*(base.md))))...
      -X2*(e1*(base.BW2*(base.WM*(base.md)))) ...
      -X3*(e1*(base.BW3*(base.WM*(base.md)))) ...
      -X4*(e1*(base.BW4*(base.WM*(base.md)))) ...
      -X5*(e1*(base.BW5*(base.WM*(base.md)))) ...
      -X6*(e1*(base.BW6*(base.WM*(base.md)))) ...
      -X7*(e1*(base.BW7*(base.WM*(base.md)))) ...
      -X8*(e1*(base.BW8*(base.WM*(base.md))));

bak12=-X1*(e2*(base.BW1*(base.WM*(base.md)))) ...
      -X2*(e2*(base.BW2*(base.WM*(base.md)))) ...
      -X3*(e2*(base.BW3*(base.WM*(base.md)))) ...
      -X4*(e2*(base.BW4*(base.WM*(base.md)))) ...
      -X5*(e2*(base.BW5*(base.WM*(base.md)))) ...
      -X6*(e2*(base.BW6*(base.WM*(base.md)))) ...
      -X7*(e2*(base.BW7*(base.WM*(base.md)))) ...
      -X8*(e2*(base.BW8*(base.WM*(base.md))));

%%%%%%%%%%%%%%%%%% FD %%%%%%%%%%%%%%%%%%%%%%%%%

X1=(sqrt(-1)*base.AM1*f+base.AD1+sqrt(-1)*base.AL1/f);
X2=(sqrt(-1)*base.AM2*f+base.AD2+sqrt(-1)*base.AL2/f);
X3=(sqrt(-1)*base.AM3*f+base.AD3+sqrt(-1)*base.AL3/f);
X4=(sqrt(-1)*base.AM4*f+base.AD4+sqrt(-1)*base.AL4/f);
X5=(base.AD5+sqrt(-1)*base.AL5/f);
X6=(base.AD6+sqrt(-1)*base.AL6/f);
X7=(base.AD7+sqrt(-1)*base.AL7/f);
X8=(base.AD8+sqrt(-1)*base.AL8/f);

bak21=-X1*(e1*(base.AW1*(base.WM*(base.md))))...
      -X2*(e1*(base.AW2*(base.WM*(base.md)))) ...
      -X3*(e1*(base.AW3*(base.WM*(base.md)))) ...
      -X4*(e1*(base.AW4*(base.WM*(base.md)))) ...
      -X5*(e1*(base.AW5*(base.WM*(base.md)))) ...
      -X6*(e1*(base.AW6*(base.WM*(base.md)))) ...
      -X7*(e1*(base.AW7*(base.WM*(base.md)))) ...
      -X8*(e1*(base.AW8*(base.WM*(base.md))));

bak22=-X1*(e2*(base.AW1*(base.WM*(base.md)))) ...
      -X2*(e2*(base.AW2*(base.WM*(base.md)))) ...
      -X3*(e2*(base.AW3*(base.WM*(base.md)))) ...
      -X4*(e2*(base.AW4*(base.WM*(base.md)))) ...
      -X5*(e2*(base.AW5*(base.WM*(base.md)))) ...
      -X6*(e2*(base.AW6*(base.WM*(base.md)))) ...
      -X7*(e2*(base.AW7*(base.WM*(base.md)))) ...
      -X8*(e2*(base.AW8*(base.WM*(base.md))));

bak1=bak11+bak21;
bak2=bak12+bak22;

else
error('P*vec hata var\n');

end

end

