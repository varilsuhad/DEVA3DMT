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






% norm(matrv1)


% min(matrc)
% max(matrc)
% 
% v1=csr2sparse(matrv1,matrr, matrc,pas.N);
% v2=csr2sparse(matrv2,matrr, matrc,pas.N);
% 
% WL=csr2sparse(pas.WLval{ii},pas.WLrow{ii}, pas.WLcol{ii},pas.WLC{ii});
% 
% bak1=v1.'*(WL'*dd1);
% bak2=v2.'*(WL'*dd2);


[bak1,bak2]=Lmergemex2P(pas.WLrow{ii}, pas.WLcol{ii}, complex(pas.WLval{ii}), int32(pas.WLC{ii}),int32(pas.N), int32(okT), gpuArray(complex(dd1)),gpuArray(complex(dd2)),matrr,matrc,matrv1,matrv2);


% WLr=pas.WLrow{1};
% WLc=pas.WLcol{1};
% WLr=pas.WLval{1};
% WLC=pas.WLC{1};
% WRr=pas.WRrow{1};
% WRc=pas.WRcol{1};
% WRv=pas.WRval{1});
% WRC=pas.WRC{1};
% WRsize=pas.WRsize;
% C1r=pas.C1row;
% C1c=pas.C1col;
% C1v=(pas.C1val);
% C2r=pas.C2row;
% C2c=pas.C2col;
% C2v=(pas.C2val);
% C3r=pas.C3row;
% C3c=pas.C3col;
% C3v=(pas.C3val);
% C4r=pas.C4row;
% C4c=pas.C4col;
% C4v=(pas.C4val);
% C5r=pas.C5row;
% C5c=pas.C5col;
% C5v=(pas.C5val);
% C6r=pas.C6row;
% C6c=pas.C6col;
% C6v=(pas.C6val);
% 
% save('Cmergemexdeb.mat','Z0','PA0',


if(pas.WRC{ii}>0)
    [bak3]=Cmergemex(Z0, complex(PA0), pas.WLrow{ii}, pas.WLcol{ii}, complex(pas.WLval{ii}), int32(pas.WLC{ii}),pas.WRrow{ii}, pas.WRcol{ii}, ...
        complex(pas.WRval{ii}), int32(pas.WRC{ii}), int32(pas.WRsize), int32(okT), ...
        pas.C1row,pas.C1col,complex(pas.C1val), int32(pas.CC), pas.C2row, pas.C2col, complex(pas.C2val), pas.C3row, pas.C3col, complex(pas.C3val), ... 
        pas.C4row, pas.C4col, complex(pas.C4val), pas.C5row, pas.C5col, complex(pas.C5val), pas.C6row, pas.C6col, complex(pas.C6val),gpuArray(complex(dd3)), ...
        int32(ii),int32(freqdep),int32(complexdist));
% bak3=0;

    %     [C] = Cform(pas,veri,set,ii,length(pas.flist));
    % %     bak30;
    %     if( okT==1)
    %     bak3=C*dd3;
    %     elseif(okT==2)
    % %         size(C)
    % %         size(dd3)
    %     bak3=transpose(C)*dd3;
    %     % bak3=(C')*dd3;
    %     end
else
    if(okT==1)
    bak3=gpuArray(zeros(size(bak1)));
    else
    bak3=zeros(0,1);
    end
end


end