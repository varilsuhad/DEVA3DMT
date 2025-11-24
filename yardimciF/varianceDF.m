% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
% Estimate variances for impedance-derived determinant data and handle NaN/zero entries.
function [det,sigD] = varianceDF(Zd,e)

    al=sum(isnan(Zd(:)));

    if(al>0)
        det=NaN(1)+sqrt(-1)*NaN(1);
        sigD=det;
        return
    end

det=sqrt((Zd(1,1)^2+Zd(1,2)^2+Zd(2,1)^2+Zd(2,2)^2)/2);

zxx1=0.5*inv(det)*Zd(1,1);
zxy1=0.5*inv(det)*Zd(1,2);
zyx1=0.5*inv(det)*Zd(2,1);
zyy1=0.5*inv(det)*Zd(2,2);

sigD(1,1)=abs(sqrt(zxx1^2*e(1,1)^2+zxy1^2*e(1,2)^2 +zyx1^2*e(2,1)^2+zyy1^2*e(2,1)^2));

end

