% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [dz]=anablok3DMDF(f,set)


oran=set.minblok;
ks=set.dzks;

rort=set.sdro;
d=1.5*503.2*sqrt(rort./f)';

mi=min(503.2*sqrt(rort./f)');
ma=max(d);

if(set.autosd==0)
    
ma=set.manualsd;   
end


st=mi/oran; %80 iyi gibi


c=1;
dz(1)=st;
tot=dz(1);
while(tot<ma)
c=c+1;
dz(c)=dz(c-1)*ks;
tot=tot+dz(c);
end

% ro=ones(1,1,length(dz))*rort;
dz=dz';
fprintf('Maximum Depth= %.1fkm and minimum block height(sd) = %.2fm\n z-dir block number=%d\n',ma/1000,st,c);

end