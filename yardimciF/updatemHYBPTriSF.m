% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [base] = updatemHYBPTriSF(base,set,dm)


ara=log(base.m)+dm(1:base.totP);
base.m=exp(ara);
nf=length(base.WE.WR);

if (base.totC>0)
ara=dm(base.totP+1:end);
base.D1=base.D1+ara;

if(set.complexdist==1)
ara=ara(1:end/2)+sqrt(-1)*ara(end/2+1:end);
end
[r1,r2]=size(base.D{1});  
if (set.freqdependent==0)

    for i=1:nf  
    ekle=base.WE.WR{i}*ara;    
    ekle=transpose(reshape(ekle,[r2 r1]));    
    base.D{i}=base.D{i}+ekle;   
    end
else 
    son=0;    
    for i=1:nf  
    ek=size(base.WE.WR{i},2);    
    ekle=base.WE.WR{i}*ara(son+1:son+ek);    
    ekle=transpose(reshape(ekle,[r2 r1]));    
    base.D{i}=base.D{i}+ekle;   
    son=son+ek;
    end     
end


end



end

