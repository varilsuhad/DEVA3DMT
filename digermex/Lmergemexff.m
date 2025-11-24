% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [bak1,bak2,bak3] = Lmergemexff(pas,veri,ii,okT,dd1,dd2,dd3,set)

Z0=gpuArray(veri.Z0);
D=gpuArray(pas.D{ii});
T0=gpuArray(veri.T0);
A0=gpuArray(veri.A0);
IH=gpuArray(veri.IH);
IE=gpuArray(veri.IE);
DPT=gpuArray(veri.DPT);
PT=gpuArray(veri.PT);
PV=gpuArray(veri.PV);
TA=gpuArray(veri.TA);
T=gpuArray(veri.T);
DPV=gpuArray(veri.DPV);
DPB=gpuArray(veri.DPB);
PTV=gpuArray(veri.PTV);
PA0=gpuArray(veri.PA0);
Z=gpuArray(veri.Z);
DD=gpuArray(veri.Det);

f=veri.f;

freqdep=set.freqdependent;
complexdist=set.complexdist;

[matrr,matrc,matrv1,matrv2]=Lmergemex1P(double(f), complex(D), Z0, T0, T, A0, IH, IE, complex(DPT), complex(PT), complex(PV), ...
    complex(TA), complex(DPV), complex(DPB), complex(PTV), complex(PA0), complex(Z),complex(DD),...
    pas.L1xrow,pas.L1xcol,pas.L1xval,pas.L1yrow,pas.L1ycol,pas.L1yval,pas.L2xrow,pas.L2xcol,complex(pas.L2xval), ...
    pas.L2yrow, pas.L2ycol, complex(pas.L2yval), pas.Hxrow, pas.Hxcol, complex(pas.Hxval), pas.Hyrow, pas.Hycol, complex(pas.Hyval), ...
pas.Hzrow, pas.Hzcol, complex(pas.Hzval),int32(pas.N));

[bak1,bak2]=Lmergemex2P(pas.WLrow{ii}, pas.WLcol{ii}, complex(pas.WLval{ii}), int32(pas.WLC{ii}),int32(pas.N), int32(okT), gpuArray(complex(dd1)),gpuArray(complex(dd2)),matrr,matrc,matrv1,matrv2);

if(pas.WRC{ii}>0)
[bak3]=Cmergemex(Z0, complex(PA0), pas.WLrow{ii}, pas.WLcol{ii}, complex(pas.WLval{ii}), int32(pas.WLC{ii}),pas.WRrow{ii}, pas.WRcol{ii}, ...
        complex(pas.WRval{ii}), int32(pas.WRC{ii}), int32(pas.WRsize), int32(okT), ...
        pas.C1row,pas.C1col,complex(pas.C1val), int32(pas.CC), pas.C2row, pas.C2col, complex(pas.C2val), pas.C3row, pas.C3col, complex(pas.C3val), ...
        pas.C4row, pas.C4col, complex(pas.C4val), pas.C5row, pas.C5col, complex(pas.C5val), pas.C6row, pas.C6col, complex(pas.C6val),gpuArray(complex(dd3)), ...
        int32(ii),int32(freqdep),int32(complexdist));
else
    if(okT==1)
    bak3=gpuArray(zeros(size(bak1)));
    else
    bak3=zeros(0,1);
    end
end

end
