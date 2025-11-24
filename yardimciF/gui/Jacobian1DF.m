% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [J] = Jacobian1DF(ro,dz,f,set)

[Z,go1,faz1] = MT1DF(dz,ro,f);

if (set.MT1DInvEmp==0)

    hep1=[log(go1);faz1];
    for i=1:length(ro)
        ro2=ro;
        s1=ro2(i);
        s2=ro2(i)*1.0001;
        h=log(s2)-log(s1);
        ro2(i)=s2;

        [~,go2,faz2] = MT1DF(dz,ro2,f);
        hep2=[log(go2);faz2];
        J(:,i)=[hep2-hep1]/h;
    end

elseif(set.MT1DInvEmp==1)
    hep1=Z;
    for i=1:length(ro)
        ro2=ro;
        s1=ro2(i);
        s2=ro2(i)*1.0001;
        h=log(s2)-log(s1);
        ro2(i)=s2;

        [Z2,go2,faz2] = MT1DF(dz,ro2,f);
        hep2=Z2;
        J(:,i)=[hep2-hep1]/h;
    end       
else
    
    fprintf('Bir hata var 1D inversion \n');
end


end

