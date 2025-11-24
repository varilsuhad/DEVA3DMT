% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [grad,set,PC] = JTdd3DMTHYBMDTriGSF(PC,set)

aa=tic;

top=0;
top2=0;

grad=zeros(PC.totP+PC.totC,1,'gpuArray');

nf=length(PC.WE.Wd);
res=zeros(nf,2);
for i=1:nf

    f=PC.veri(i).f;

    if (set.frekans==1)
        fprintf('Frequency=%f \n',f);
    end

    [PC] = freqmergeCF(PC,PC.veri(i).f,set);

    A1=(PC.WE.Wd{i});
    A2=(PC.veri(i).ddi);
    dd=A1*(A1*conj(A2));
    PC.veri(i).wdd=dd;
    dd=gpuArray(dd);
    [bb1,bb2,bb3] = Lmergemexff(PC,PC.veri(i),i,2,dd,dd,dd,set);


    if(f>=10)
        ek=10^-1;
    else
        ek=1;
    end

    if(set.zeroinitialx==1)
        PC.xJT{i}=zeros(size(PC.xJT{i}));
    end

    ss=tic;
    b12=[bb1;bb2];

    if(set.SPpreconditioning==1)
        if(set.MatrixStacking==1)
            [xx,r,relres]=gpbicgmexP2FFsi(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,(PC.valM),set.limDiger*ek,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        else
            [xx,r,relres]=gpbicgmexP2FFsS(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,(PC.valM),set.limDiger*ek,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        end
    else
        if(set.MatrixStacking==1)
            [xx,r,relres]=gpbicgmexP2FFdi(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,(PC.valM),set.limDiger*ek,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        else
            [xx,r,relres]=gpbicgmexP2FFsD(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,(PC.valM),set.limDiger*ek,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        end
    end

    if(set.relres==1)
        fprintf('relres=%e iter=%d\n',relres,length(r));
    end
    if(relres>set.limDiger*ek)
        fprintf('warning -> relres=%e but the limit is set to %e, frequency=%fHz\n',relres,set.limDiger*ek,f);
    end

    if(relres>set.limDiger*ek)
        fprintf('solving again due to very low accuracy\n');
        PC.xJT{i}=gpuArray(zeros(size(PC.xJT{i})));
        if(set.MatrixStacking==1)
            [xx2,r2,relres2]=gpbicgmexP2FFdi(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,PC.valM,set.limForward,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        else
            [xx2,r2,relres2]=gpbicgmexP2FFsD(PC.rowA,PC.colA,PC.valA,complex(b12),PC.rowM,PC.colM,PC.valM,set.limForward,set.maxit,set.stagdetect,complex(PC.xJT{i}));
        end
        if(relres2<relres)
            fprintf('New relres=%e and the limit=%e, f=%f\n',gather(relres2),set.limDiger*ek,f);
            xx=xx2;
            relres=relres2;
            r=r2;
        else
            fprintf('Using old relres=%e and the limit=%e, f=%f\n',gather(relres),set.limDiger*ek,f);
        end
    end


    PC.rowA=[];PC.colA=[];PC.valA=[];PC.rowM=[];PC.colM=[];PC.valM=[];


    PC.xJT{i}=xx;
    xx=reshape(xx,[],2);
    x1=xx(:,1);
    x2=xx(:,2);
    PC.veri(i).s1a=(x1);   %% ikinci türevler için
    PC.veri(i).s1b=(x2);

    res(i,1:2)=relres;
    top2=top2+length(r);


    top=top+toc(ss);
    [ Px1,Px2 ] = Pmerge3DMTHYBS3MULTFF3mex(PC,PC.veri(i),2,x1,x2,set);



    ara1=(real(Px1)+real(Px2));
    grad=grad+[ara1;real(bb3)];


end
grad=gather(grad);

aa=toc(aa);
set.JTddSaveTime(end+1)=aa;

if(isfield(set,'relresJTdd')==0)
    set.relresJTdd={};
end
if(isfield(set,'JTiteration')==0)
    set.JTiteration=[];
end

set.relresJTdd{end+1}=res;
set.JTiteration(end+1)=top2;

if(set.sureJTdd==1)
    fprintf('JT*deltad is calculated in %.2f secs \n but Iteration time is %.2f secs\n Total iteration = %d\n',aa,top,top2);
end



end
