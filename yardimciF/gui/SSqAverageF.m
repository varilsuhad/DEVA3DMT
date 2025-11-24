% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [Zssq,Vssq] = SSqAverageF(data)

[nf,ni,~]=size(data);

for i=1:nf
    c=0;
    clear Z1 V1

    for j=1:ni

    Z=data(i,j,1:4);

    zxx=isnan(Z(1));
    zxy=isnan(Z(2));
    zyx=isnan(Z(3));
    zyy=isnan(Z(4));

    bol=2;

%     if(zxx==1 || zxy==1 || zyx==1 || zyy==1)
%      continue;
%     else
%      v=[data(i,j,7);data(i,j,8);data(i,j,9);data(i,j,10)];
%     end

    v=[];
    if (zxy==1)
        bol=bol-1;
%         Z(1)=NaN;
    else
        v1=data(i,j,8);
        if(~(isnan(v1) || v1==0))
        v=[v;v1];
        end
    end

    if (zyx==1)
        bol=bol-1;
%         Z(4)=NaN;
    else
        v1=data(i,j,9);
        if(~(isnan(v1) || v1==0))
        v=[v;v1];
        end
    end

    if (zxx==0)
        v1=data(i,j,7);
        if(~(isnan(v1) || v1==0))
        v=[v;v1];
        end
    end

    if (zyy==0)
        v1=data(i,j,10);
        if(~(isnan(v1) || v1==0))
        v=[v;v1];
        end
    end

    if(bol==0)
        continue;
    end
    Z(isnan(Z))=0;

    c=c+1;
    Z1(c,1)=sqrt((Z(1)^2+Z(2)^2+Z(3)^2+Z(4)^2)/bol);
    V1(c,1)=sqrt(transpose(v)*v/length(v));
    end

    Zssq(i,1)=GeometricAverageF(Z1);
    Vssq(i,1)=GeometricAverageF(V1);
end

Zssq=Zssq(:);
Vssq=Vssq(:);

end

