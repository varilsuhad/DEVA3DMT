function [base,F,set,WF,ddi,wddi] = deltad3DMTMDSF(base,set,uu)

ddi=0;
wddi=0;
nf=length(base.WE.Wd);

ind(1)=0;
for i=1:nf
r1=size(base.WE.Wd{i},1);
ind(i+1)=ind(i)+r1;    
end
F=zeros(ind(i+1),1);   
WF=zeros(ind(i+1),1);        

top=0;
for i=1:nf      
    tot=0;
    for ii=1:i-1
    tot=tot+size(base.WE.Wd{ii},1);    
    end
    ek=size(base.WE.Wd{i},1);        
    base.veri(i).ddi=base.d(tot+1:tot+ek)-base.veri(i).Fi;  
    ara=base.WE.Wd{i}*base.veri(i).ddi;
    top=top+ara'*ara;
    F(ind(i)+1:ind(i+1))=base.veri(i).Fi;
    WF(ind(i)+1:ind(i+1))=ara;        
end  


dd=nnz(imag(F))+nnz(real(F));
set.RMSall(uu)=real(sqrt(top/dd));
set.MISFIT(uu)=real(top);


end

