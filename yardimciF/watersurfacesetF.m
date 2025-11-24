% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [NK,ro] = watersurfacesetF(NK,ebhava,z,ro,set)


[ny,nx,nz,~]=size(NK);

ny1=ny-1;
nx1=nx-1;

[ny2,nx2,nz2]=size(ro);

ekbloky=(ny1-ny2)/2;
ekblokx=(nx1-nx2)/2;


% al=squeeze(NK(83,1,ebhava+1:end,3));
% goal=z(83,1);



% ro2=zeros(ny-1,nx-1,nz-1);

for i=1:nx
    for j=1:ny
    goal=z(j,i);

    if(goal>=0)
        continue;
    end

    al=squeeze(NK(j,i,ebhava+1:end,3));
    if(al(1)~=0)
        if(abs(al(1))<10^-5)
            al(1)=0;
        else
        error('watersurfacesetF surface not 0');
        end
    end

    [al,ind] = watersubF(al,goal);
    NK(j,i,ebhava+1:end,3)=al;
    if( i>=ekblokx+1 && i<ekblokx+nx2 && j>=ekbloky+1 && j<=ekbloky+ny2)
    ro(j-ekbloky,i-ekblokx,1:ind-1)=-2;
    end

    end
end



















% [ny,nx]=size(z);
% 
% c=0;c2=0;
% tol=0.00001;
% if(isfield(set,'denizblok')==0)
% set.denizblok=0;
% end
% 
% denizblok=set.denizblok;
% 
% sakla=[];
% for i=1:nx
%     for j=1:ny
%         al1=NK(j,i,ebhava+1,3); 
%         al2=z(j,i);
%         if((al1+al2)>tol)
%             error('Topology addition error\n');
%         end
% 
%         if(al1>tol)
%             c=c+1;
%             al=squeeze(NK(j,i,1:ebhava+1,3)); 
%             [geri,db]=suayarF(al,denizblok);  %%bura
% %             geri2=suayar2F(al,denizblok);  %%bura
%             NK(j,i,1:ebhava+1,3)=geri;
%             ind1=1:ebhava+1-db;
%             al1=squeeze(NK(j,i,ind1,3));
% 
%             kat=1.0;
%             g1=logaritmikkaydirF(al1,kat);
%             while(abs(g1(end-1))>300)
%             kat=kat+0.05;    
%             g1=logaritmikkaydirF(al1,kat);                
%             end
%             if(kat~=1.0)
%                c2=c2+1; 
%             end
% 
%             sakla=[sakla;kat];
% 
%             NK(j,i,ind1,3)=g1;
%         end        
% 
%     end
% end
% fprintf('%d sea floor correction\n',c);
% fprintf('%d sea-air correction\n',c2);
% 
% 
% if(c>0)
%     kat=mean(sakla);   
%     for i=1:nx
%         for j=1:ny
%             al1=NK(j,i,ebhava+1,3); 
%             al2=z(j,i);
%             if((al1+al2)>tol)
%                 error('Topo addition error\n');
%             end
%             if(al1<=tol)
%                 al1=squeeze(NK(j,i,1:ebhava+1,3));
%                 g1=logaritmikkaydirF(al1,kat);        
%                 NK(j,i,1:ebhava+1,3)=g1;
%             end        
%         end
%     end    
% end






end

