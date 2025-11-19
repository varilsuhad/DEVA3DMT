function [al,ind] = watersubF(al,gl)


gl=abs(gl);

for i=1:length(al)-1
    n1=al(i);
    n2=al(i+1);
    if( gl>=n1 && gl<=n2)
        if( ((gl-n1)<(n2-gl)) && i~=1)
        al(i)=gl;
        ind=i;        
        break;
        else
        al(i+1)=gl;
        ind=i+1;
        break;
        end
    end
end



end