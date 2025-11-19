function [datay,fy,dzy,zmaxy,roiy,roa] = anablokDataF(data,f,set)



%%% frekansa baðlý olarak olarak dzleri oluþtur
%%% Eðer istenirse verilen datayý 1D fit et
set.autosd=1;
dz=anablok3DMDF(f,set);
if(set.initial1D==1)
    [Zssq,Wssq] = SSqAverageF(data);
%     set.MT1Dfigures=1;
    [roi,~,~]=invert1DF(Zssq,Wssq,dz(:),f,set); 
else
    roi=ones(size(dz))*set.iro;
end

%%% Skindepth hesabýna göre eriþilen derinlikleri hesapla
[zmax,roa]=skindepthaverageMDF(roi,dz,set,f);
fprintf('Medium maximum depth=%.1fkm and reachable max depth=%.1fkm\n',sum(dz)/1000,zmax(end)/1000);

%%% Eðer eriþilen derinlik daha kýsaysa dzleri tekrar oluþtur 
%%% ve tekrar 1D fit yap
if(sum(dz)>zmax(end))
%     fprintf('Medium maximum depth=%.1fkm and reachable max depth=%.1fkm\n  KIRPTIM\n',sum(dz),zmax(end));   
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
    fprintf('Ortam maxd=%f ve eriþilen maxd=%f\n',sum(dz),zmax(end));
end

%%% Eðer eriþilen derinlik eldeki meshin derinliðinden daha fazla ise
%%% dz sayýsýný arttýr.
if(sum(dz)<zmax(end))
c=0;
while(sum(dz)<zmax(end))
    dz(end+1)=dz(end)*set.dzks;
    c=c+1;
end
 roi(end+1:end+c)=roi(end);   
 fprintf('Yeni ortam maxd=%f ve eklenen blok=%d\n',sum(dz),c); 
end

%%% Eðer belirlenen maximum derinlik frekansýn eriþebileceði derinlikten
%%% kýsaysa bazý frekanslarý at

if(zmax(end)>set.maxsd)
    fprintf('eriþilen maxd=%f ve berirlenen max=%f d\n',zmax(end),set.maxsd);
       
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
    fprintf('Atýlan frekans=%d ve dz blok=%d adet\n',length(ind),length(ind2));
else
    datay=data;    
    dzy=dz;    
    fy=f;
    zmaxy=zmax;
    roiy=roi;
end


end

