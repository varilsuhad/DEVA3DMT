% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [zmax,roa,roi] = skindepthaverageMDF(roi,dz,set,f)

minf=min(f);
sig=1./roi;
sig=sig(:);
dz=dz(:);

for i=1:length(roi)
asig(i)=sum(sig(1:i).*dz(1:i))/sum(dz(1:i));
end

for j=1:length(f)
minf=f(j);

% for i=1:length(roi)
%     ara1(i)=503.3*1.5*sqrt(1/asig(i)/minf);
%     ara2(i)=sum(dz(1:i));
%     if(ara2(i)>ara1(i))
%         break
%     end
% end

c=0;
son=length(asig);
while(1)
    c=c+1;
    if (c<=son)
        ara1(c)=503.3*1.5*sqrt(1/asig(c)/minf);
        ara2(c)=sum(dz(1:c));
    else
        ara1(c)=503.3*1.5*sqrt(1/asig(end)/minf);
        ara2(c)=sum(dz(1:end))+sum(dz(end)*set.dzks.^(1:(c-son)));
        roi(c)=roi(end);
    end

    if(ara2(c)>ara1(c))
        break
    end
end

if (c<=son)
    uu=min(c+2,length(dz));
    roa(j,1)=1/asig(c);
    zmax(j,1)=sum(dz(1:uu));
else
    roa(j,1)=1/asig(end);
    zmax(j,1)=503.3*1.5*sqrt(1/asig(end)/minf);
end

% zmax(j,1)=sum(dz(1:uu));
% roa(j,1)=1/asig(i);

%
% zmax(j,1)=sum(dz(1:uu));
% roa(j,1)=1/asig(i);

end

end

