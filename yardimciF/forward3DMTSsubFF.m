% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [veri,base] = forward3DMTSsubFF(base,f,set,nof)

    if(set.zeroinitialx==1)
    base.xfor{nof}=zeros(size(base.xfor{nof}));
    end

    if(f>=10)
    ek=10^-1;
    else
    ek=1;
    end

    if(set.SPpreconditioning==1)
        if(set.MatrixStacking==1)
        [x,r,relres]=gpbicgmexP2FFsi(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,(base.valM),set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        else
        [x,r,relres]=gpbicgmexP2FFsS(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,(base.valM),set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        end
    else
        if(set.MatrixStacking==1)
        [x,r,relres]=gpbicgmexP2FFdi(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,(base.valM),set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        else
        [x,r,relres]=gpbicgmexP2FFsD(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,(base.valM),set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        end
    end

    if(relres>set.limForward*ek)
    fprintf('warning -> relres=%e but the limit is set to %e, frequency=%fHz\n',relres,set.limForward*ek,f);
    end

    if(relres>set.limForward*ek)
        fprintf('solving again due to very low accuracy\n');
        base.xfor{nof}=gpuArray(zeros(size(base.xfor{nof})));
        if(set.MatrixStacking==1)
        [x2,r2,relres2]=gpbicgmexP2FFdi(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,base.valM,set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        else
        [x2,r2,relres2]=gpbicgmexP2FFsD(base.rowA,base.colA,base.valA,complex(base.bstack),base.rowM,base.colM,(base.valM),set.limForward*ek,set.maxit,set.stagdetect,complex(base.xfor{nof}));
        end
        if(relres2<relres)
        fprintf('New relres=%e and the limit=%e, f=%f\n',gather(relres2),set.limForward*ek,f);
        x=x2;
        relres=relres2;
        r=r2;
        else
        fprintf('Using old relres=%e and the limit=%e, f=%f\n',gather(relres),set.limForward*ek,f);
        end

    end

    base.rowA=[];
    base.colA=[];
    base.valA=[];
    base.rowM=[];
    base.colM=[];
    base.valM=[];

x=gather(reshape(x,[],2));
res(1,1:2)=gather(relres);
veri=verihesapF( base,x,f,nof,res,set);
veri.rr=r;

end

