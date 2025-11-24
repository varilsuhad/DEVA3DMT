% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Insert a water layer into the mesh coordinates and resistivity model according to bathymetry.
function [NK,ro] = watersurfacesetF(NK,ebhava,z,ro,set)

[ny,nx,nz,~]=size(NK);

ny1=ny-1;
nx1=nx-1;

[ny2,nx2,nz2]=size(ro);

ekbloky=(ny1-ny2)/2;
ekblokx=(nx1-nx2)/2;

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

end

