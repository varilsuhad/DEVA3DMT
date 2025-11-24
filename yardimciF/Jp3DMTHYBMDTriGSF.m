function [grad,set,PC] = Jp3DMTHYBMDTriGSF(PC,p,set)

aa=tic;


top=0;
top2=0;

nf=length(PC.WE.Wd);
totP=PC.totP;  
ind=zeros(nnz(PC.flist)+1,1);
c=0;
for i=1:nf
    c=c+1;
    r1=size(PC.WE.Wd{i},1);
    ind(c+1)=ind(c)+r1;    
end
grad=(zeros(ind(c+1),1));
res=zeros(nf,2);
c=0;
for i=1:nf
       
    ok=PC.flist(i);
    if(ok==0)
        continue;
    end     
    c=c+1;
    f=PC.veri(i).f;    

    if (set.frekans==1)
    fprintf('frekans=%f \n',f);
    end

      [PC] = freqmergeCF(PC,PC.veri(i).f,set);
      [bb1,bb2] = Pmerge3DMTHYBS3MULTFF3mex(PC,PC.veri(i),1,p(1:totP),p(1:totP),set);  
        

    if(f>=10)
        ek=10^-1;
    else
        ek=1;
    end

     ss=tic;         

     bb=[(bb1);(bb2)];

     PC.veri(i).bb1=bb1;
     PC.veri(i).bb2=bb2;
     

    if(set.zeroinitialx==1)
    PC.xJp{i}=zeros(size(PC.xJp{i}));
    end     

    if(set.SPpreconditioning==1)
        if(set.MatrixStacking==1)                              
         [xx,r,relres]=gpbicgmexP2FFsi(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,(PC.valM),set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));   
        else
         [xx,r,relres]=gpbicgmexP2FFsS(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,(PC.valM),set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));   
        end
    else
        if(set.MatrixStacking==1)                                      
         [xx,r,relres]=gpbicgmexP2FFdi(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,(PC.valM),set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));        
        else
         [xx,r,relres]=gpbicgmexP2FFsD(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,(PC.valM),set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));        
        end
    end

    if(set.relres==1)        
    fprintf('relres=%e iter=%d\n',relres,length(r));
    end
    if(relres>set.limDigerGN*ek)
    fprintf('warning -> relres=%e but the limit is set to %e, frequency=%fHz\n',relres,set.limDigerGN*ek,f);    
    end   

    % if(relres>10^-8)
    if(relres>set.limDigerGN*ek)

    fprintf('solving again due to very low accuracy\n');
        PC.xJp{i}=gpuArray(zeros(size(PC.xJp{i})));
        if(set.MatrixStacking==1)                                          
        [xx2,r2,relres2]=gpbicgmexP2FFdi(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,PC.valM,set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));   
        else
        [xx2,r2,relres2]=gpbicgmexP2FFsD(PC.rowA,PC.colA,PC.valA,complex(bb),PC.rowM,PC.colM,PC.valM,set.limDigerGN*ek,set.maxit,set.stagdetect,complex(PC.xJp{i}));   
        end
        if(relres2<relres)
        fprintf('New relres=%e and the limit=%e, f=%f\n',gather(relres2),set.limDigerGN*ek,f);
        xx=xx2;
        relres=relres2;
        r=r2;
        else
        fprintf('Using old relres=%e and the limit=%e, f=%f\n',gather(relres),set.limDigerGN*ek,f);
        end
    end   

    PC.rowA=[];
    PC.colA=[];
    PC.valA=[];
    PC.rowM=[];
    PC.colM=[];
    PC.valM=[];     

    PC.xJp{i}=xx;          
    xx=reshape(xx,[],2);
    x1=(xx(:,1));
    x2=(xx(:,2));

    PC.veri(i).s2a=(x1);   %% ikinci türevler için  
    PC.veri(i).s2b=(x2);  

    res(i,1:2)=relres;        
    top=top+toc(ss); 
    top2=top2+length(r); 


    [Lx1,Lx2,Cp] = Lmergemexff(PC,PC.veri(i),i,1,x1,x2,p(totP+1:end),set);
    ara1=Lx1+Lx2;
    grad(ind(c)+1:ind(c+1))=gather(ara1+Cp);

    PC.veri(i).pddi=gather(ara1+Cp);

end


if(isfield(set,'relresJp')==0)
    set.relresJp={};
end

if(isfield(set,'JpSaveTime')==0)
    set.JpSaveTime=[];
end

if(isfield(set,'Jpiteration')==0)
    set.Jpiteration=[];
end


aa=toc(aa);
set.JpSaveTime(end+1)=aa;
set.relresJp{end+1}=res;
set.Jpiteration(end+1)=top2;


if(set.sureJp==1)
fprintf('J*p is calculated in %.2f secs\n but Iteration time is %.2f secs\n Total iteration = %d \n',aa,top,top2);
end



end

