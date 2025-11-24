% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Persist intermediate inversion outputs (fields, residuals, models) to disk for later inspection.
function [] = saveVars3DMTF(base,set,F,ii,path0,ok,r0,s1,y1,ro1)

if (ok==1)
bas.m=base.m;
bas.D=base.D;
bas.f=base.f;
bas.ekblokx=base.ekblokx;
bas.ekbloky=base.ekbloky;
bas.ebhava=base.ebhava;
bas.recv=base.recv;
bas.EL=base.EL;
bas.NK=base.NK;
bas.brecvlist=base.brecvlist;
bas.d=base.d;
bas.smooth=base.smooth;
bas.veriname=base.veriname;
bas.meshname=base.meshname;
bas.Adx=base.Adx;
bas.Ady=base.Ady;
bas.liste=base.liste;
bas.totC=base.totC;
bas.totD=base.totD;
bas.totP=base.totP;
bas.WE.WL=base.WE.WL;
bas.WE.WR=base.WE.WR;
bas.WE.WD=base.WE.WD;
bas.ro=base.ro;
bas.FE=base.FE;
bas.FD=base.FD;
bas.yuzey=base.yuzey;

if(isfield(base,'stationid'))
bas.stationid=base.stationid;
end

if(isfield(base,'koor'))
bas.koor=base.koor;
bas.koorm=base.koorm;
end

base=bas;
dir =strcat(path0,'.mat');
save(dir,'base','set');
% save(dir,'bas','set');

end

if(ok==2)
    dir =strcat(path0,num2str(ii),'.mat');
    D=base.D;
    D1=base.D1;
    m=base.m;

    % save(dir,'m','F','D','D1','set');

    smooth=base.smooth;
    save(dir,'m','F','D','D1','set','smooth');

end

if(ok==3)
    dir =strcat(path0,num2str(0),'.mat');
    ML=base.ML;
    MU=base.MU;
    save(dir,'r0','s1','y1','ro1','ML','MU');
end

end

