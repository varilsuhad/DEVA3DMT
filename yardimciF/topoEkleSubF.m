function [aa] = topoEkleSubF(nhz,goal,ebhava,ek)
%UNTÝTLED Summary of this function goes here

% ek=0.5;

kati=1+ek;
katg=1-ek;

if(goal==0)
aa=[0 cumsum(nhz)']';
ok=0;
elseif(goal>0)
ok=1;
top=0;    
nhzy=nhz;
for i=ebhava:-1:1
    
    if( i==2)
        error('topografya hatasi1');
    end    
    
    top=top+nhzy(i)*(1-katg);
    nhzy(i)=nhzy(i)*katg;
    if (top>goal)
    fark=top-goal;
    nhzy(i)=nhzy(i)+fark;    
    break;
    end

end

top=0;
for i=ebhava+1:1:length(nhz)
    
    if( i==length(nhz)-1)
        error('topografya hatasi2');
    end      
    
    top=top+nhzy(i)*(kati-1);
    nhzy(i)=nhzy(i)*kati;
    if (top>goal)
    fark=top-goal;
    nhzy(i)=nhzy(i)-fark;    
    break;
    end
  
end
aa=[0 cumsum(nhzy)']';
elseif (goal<0)
ok=1;    
top=0;    
nhzy=nhz;
for i=ebhava:-1:1
    
    if( i==2)
        error('topografya hatasi3');
    end    
    
    top=top+nhzy(i)*(kati-1);
    nhzy(i)=nhzy(i)*kati;
    if (top>abs(goal))
    fark=top-abs(goal);
    nhzy(i)=nhzy(i)-fark;    
    break;
    end

end

top=0;
for i=ebhava+1:1:length(nhz)
    
    if( i==length(nhz)-1)
        error('topografya hatasi4');
    end       
    
    top=top+nhzy(i)*(1-katg);
    nhzy(i)=nhzy(i)*katg;
    if (top>abs(goal))
    fark=top-abs(goal);
    nhzy(i)=nhzy(i)+fark;    
    break;
    end
 
end    
aa=[0 cumsum(nhzy)']';
else
    error('Topografya bir hata\n');
end


end

