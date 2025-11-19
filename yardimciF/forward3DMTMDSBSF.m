function [PC,set,veri] = forward3DMTMDSBSF(PC,f,set)

aa=tic;


top=0;
top2=0;
for i=1:length(f)
    [PC] = freqmergeCF(PC,f(i),set);
    ss=tic;
    [veri(i),PC]=forward3DMTSsubFF(PC,f(i),set,i);
    top=top+toc(ss);
    top2=top2+length(veri(i).rr);
    PC.xfor{i}=[veri(i).e1;veri(i).e2];
end    
PC.veri=veri;


aa=toc(aa);
set.forwardSaveTime(end+1)=aa;

if(set.sureForward==1)
fprintf('Forward modeling with finer mesh is completed in %.2f secs \n but Iteration time is %.2f secs\n Total iteration no = %d \n\n',aa,top,top2);
end


res=zeros(length(f),2);   
for i=1:length(f)
    res(i,:)=veri(i).res;
end       


if(isfield(set,'relresForward')==0)
    set.relresForward={};
end
if(isfield(set,'forwarditerationA')==0)
    set.forwarditerationA=[];
end

set.relresForward{end+1}=res;
set.forwarditerationA(end+1)=top2;




end

