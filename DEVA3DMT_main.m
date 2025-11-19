function [] = DEVA3DMT_main(inputmatrix,outputfolder,outputname)
format longG

%%%%% Ayar Yükle
[set]=settings3DMTF();


d=gpuDevice;
reset(d)
d.CachePolicy='minimum';


load(inputmatrix);
set.olcurotasyonxy=aci;
fprintf("\nRotation angle=%.2f degrees\n",aci);

[base,set] = basekurSF(x,y,z,f,recv,set,data);  



% % 100 ohm.m starting model
% base.m=ones(length(base.m),1)*(1/100);


pause(0.1);
path=strcat(outputfolder,'\',outputname);

if(exist('stationid'))
base.stationid=stationid;
end
if(exist('koor'))
base.koor=koor;
base.koorm=koorm;
end

if(exist('Adx')==1)
base.Adx=Adx;
base.Ady=Ady;
end

if(exist('veriname')==0)
veriname='veri';
end
if(exist('meshname')==0)
meshname='meshname';
end

base.veriname=veriname;
base.meshname=meshname;



saveVars3DMTF(base,set,[],[],path,1);

fprintf("\n Starting...\n");
start=tic;

[base] = basemergeMexSF(base,set);
[base,set] = forward3DMTMDSBSF(base,f,set);

[base,F,set] = deltad3DMTMDSF(base,set,1);   
fprintf('\nRMS=%e misfit=%e \n',set.RMSall(1),set.MISFIT(1));


for uu=1:set.maxitInv
aa=tic;
set.rs=set.rs+1;

[JTdd,set,base]=JTdd3DMTHYBMDTriGSF(base,set);
r=2*JTdd-2*base.CC*[log(base.m);(base.D1-base.D0)]; %reel r

m=50;
r0=r;
e(1,1)=norm(r0);
v(:,1)=r0/e(1,1);

MA=base.M;

U=ichol(MA);
base.ML=transpose(U);
base.MU=U;
clear sav

for j=1:m

if(set.zeroinitialx==0 && j==1)
set.zeroinitialx=1;
modded=1;
else
modded=0;
end

    ara1=base.MU\(base.ML\v(:,j));
    [~,set,base]=Jp3DMTHYBMDTriGSF(base,gpuArray(ara1),set);
    [JTdd1,set,base]=JTdd3DMTHYBMDTriGSFF(base,set);

    w(:,j+1)=2*JTdd1+2*(base.CC)*ara1;    
    for i=1:j
    h(i,j)=w(:,j+1)'*v(:,i);
    w(:,j+1)=w(:,j+1)-h(i,j)*v(:,i);
    end
    
    h(j+1,j)=norm(w(:,j+1));
    v(:,j+1)=w(:,j+1)/h(j+1,j);
    
    for i=1:j-1
    ara=[c(i) s(i);-s(i) c(i)]*[h(i,j);h(i+1,j)];
    h(i,j)=ara(1);
    h(i+1,j)=ara(2);
    end
    
    c(j)=h(j,j)/sqrt(h(j,j)^2+h(j+1,j)^2);
    s(j)=h(j+1,j)/sqrt(h(j,j)^2+h(j+1,j)^2);
    e(j+1,1)=-s(j)*e(j,1);
    e(j,1)=c(j)*e(j,1);
    h(j,j)=sqrt(h(j,j)^2+h(j+1,j)^2);
    h(j+1,j)=0;
    
    ne=norm(e(j+1,1))/norm(r);
    lst(j)=ne;
    
    y=h(1:j,1:j)\e(1:j);
    dm=base.MU\(base.ML\(v(:,1:j)*y));
    
    fprintf('\n\n GN step=%d ve norm(e)=%e norm(dm)=%e \n\n',j,ne,norm(dm));

    sav(j,1)=norm(dm);
    sav(j,2)=ne;    
    sav(j,3)=norm(2*JTdd1);
    sav(j,4)=norm(2*(base.CC)*ara1);

    if(set.zeroinitialx==1 && j==1 && modded==1)
    set.zeroinitialx=0;
    modded=0;
    end

    if(ne<set.GNlimit)
        break;
    end
    % if(norm(dm)>20)
    % break;
    % end
end
clear w h lst c s e v r0 y


set.inneriterGN{uu}=sav;

alfa=1;

[base]=updatemHYBPTriSF(base,set,dm);
[base] = basemergeMexSF(base,set);
[base,set] = forward3DMTMDSBSF(base,f,set);
[base,F,set] = deltad3DMTMDSF(base,set,uu+1);   


fprintf('\nRMS=%e misfit=%e alfa=%f \n',set.RMSall(uu+1),set.MISFIT(uu+1),alfa);


[set,base] = StopOrGo3DLF(set,base,dm,alfa,uu,ne); 

aac=toc(aa);
set.OneItSaveTime(end+1)=aac;

if(set.sureOneIt==1)
fprintf('One Iteration took %.2f secs\n\n',aac); 
end

saveVars3DMTF(base,set,F,uu,path,2); 
fprintf('%s\n',datetime)

if(set.stop==1)
    fprintf('STOPPED\n');
    break
end


end

tottime=toc(start);
fprintf('\nInversion Completed in %.2f mins\n',tottime/60); 



end
