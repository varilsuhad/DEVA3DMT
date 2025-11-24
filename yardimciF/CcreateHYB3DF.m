% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Build distortion selector matrices and initialize default coefficients for all surface sites.
function [DM,D] = CcreateHYB3DF(yuzey,f)

tot=size(yuzey,1);

c=0;cc=0;

%%%Zxx(obs)
for i=1:tot
%%%Cxx
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
%%%Cxy
c=c+1;
%%%Cyx
c=c+1;
%%%Cyy
c=c+1;
%%%Czx
c=c+1;
%%%Czy
c=c+1;
end
C1=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

for i=1:tot
%%%Cxx
c=c+1;
%%%Cxy
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
%%%Cyx
c=c+1;
%%%Cyy
c=c+1;
%%%Czx
c=c+1;
%%%Czy
c=c+1;
end
C2=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

for i=1:tot
%%%Cxx
c=c+1;
%%%Cxy
c=c+1;
%%%Cyx
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
%%%Cyy
c=c+1;
%%%Czx
c=c+1;
%%%Czy
c=c+1;
end
C3=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

for i=1:tot
%%%Cxx
c=c+1;
%%%Cxy
c=c+1;
%%%Cyx
c=c+1;
%%%Cyy
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
%%%Czx
c=c+1;
%%%Czy
c=c+1;
end
C4=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

for i=1:tot
%%%Cxx
c=c+1;
%%%Cxy
c=c+1;
%%%Cyx
c=c+1;
%%%Cyy
c=c+1;
%%%Czx
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
%%%Czy
c=c+1;
end
C5=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

for i=1:tot
%%%Cxx
c=c+1;
%%%Cxy
c=c+1;
%%%Cyx
c=c+1;
%%%Cyy
c=c+1;
%%%Czx
c=c+1;
%%%Czy
c=c+1;
cc=cc+1;ix(cc)=i;iy(cc)=c;iv(cc)=1;
end
C6=sparse(ix,iy,iv,tot,tot*6);
c=0;cc=0;clear ix iy iv

DM.C1=C1;
DM.C2=C2;
DM.C3=C3;
DM.C4=C4;
DM.C5=C5;
DM.C6=C6;

% Olması gereken Distortion değerleri default values D0 için olabilir ya da
% initial
ara=[1 0 0 1 0 0];
for i=1:length(f)
D{i}(1:tot,1:6)=complex(repmat(ara,tot,1));
end

end

