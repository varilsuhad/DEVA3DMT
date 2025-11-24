% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Prepare frequency-dependent layer thicknesses and merged data arrays for the main block model.
function [datay,fy,dzy,zmaxy,roiy,roa] = anablokDataF(data,f,set)

%%% frekansa bağlı olarak olarak dzleri oluştur
%%% Eğer istenirse verilen datayı 1D fit et
set.autosd=1;
dz=anablok3DMDF(f,set);
if(set.initial1D==1)
    [Zssq,Wssq] = SSqAverageF(data);
    [roi,~,~]=invert1DF(Zssq,Wssq,dz(:),f,set);
else
    roi=ones(size(dz))*set.iro;
end

%%% Skindepth hesabına göre erişilen derinlikleri hesapla
[zmax,roa]=skindepthaverageMDF(roi,dz,set,f);
fprintf('Medium maximum depth=%.1fkm and reachable max depth=%.1fkm\n',sum(dz)/1000,zmax(end)/1000);

%%% Eğer erişilen derinlik daha kısaysa dzleri tekrar oluştur
%%% ve tekrar 1D fit yap
if(sum(dz)>zmax(end))
    set.manualsd=zmax(end);
    set.autosd=0;
    set.sdro=roa(1);
    dz=anablok3DMDF(f,set);
    if(set.initial1D==1)
        [roi,~,~]=invert1DF(Zssq,Wssq,dz(:),f,set);
    else
        roi=ones(size(dz))*set.iro;
    end
    [zmax,roa]=skindepthaverageMDF(roi,dz,set,f);
    fprintf('Ortam maxd=%f ve erişilen maxd=%f\n',sum(dz),zmax(end));
end

%%% Eğer erişilen derinlik eldeki meshin derinliğinden daha fazla ise
%%% dz sayısını arttır.
if(sum(dz)<zmax(end))
c=0;
while(sum(dz)<zmax(end))
    dz(end+1)=dz(end)*set.dzks;
    c=c+1;
end
 roi(end+1:end+c)=roi(end);
 fprintf('Yeni ortam maxd=%f ve eklenen blok=%d\n',sum(dz),c);
end

%%% Eğer belirlenen maximum derinlik frekansın erişebileceği derinlikten
%%% kısaysa bazı frekansları at

if(zmax(end)>set.maxsd)
    fprintf('erişilen maxd=%f ve belirlenen max=%f d\n',zmax(end),set.maxsd);

    fy=f;
    zmaxy=zmax;
    roiy=roi;
    dzy=dz;
    datay=data;

    ind=find(zmax>set.maxsd);
    fy(ind)=[];
    zmaxy(ind)=[];
    datay(ind,:,:)=[];

    sdz=cumsum(dz);
    ind2=find(sdz>set.maxsd);
    roiy(ind2)=[];
    dzy(ind2)=[];
    fprintf('Atılan frekans=%d ve dz blok=%d adet\n',length(ind),length(ind2));
else
    datay=data;
    dzy=dz;
    fy=f;
    zmaxy=zmax;
    roiy=roi;
end

end

