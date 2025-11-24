function [set,base] = StopOrGo3DLF(set,base,dm,alfa,uu,ne)

bas=0;

% if(uu==1)
%     set.rdmd(uu)=0;
% else
%     set.rdmd(uu)=(set.MISFIT(uu)-set.MISFIT(uu-1))/set.MISFIT(uu-1)*100;
% end


    set.rdmd(uu)=(set.MISFIT(uu+1)-set.MISFIT(uu))/set.MISFIT(uu)*100;

    set.modelr(uu)=log(base.m)'*base.C1'*base.C1*log(base.m);
    set.distr(uu)=(base.D1-base.D0)'*(base.D1-base.D0);    
    set.ndp(uu)=norm(dm(1:base.totP));
    set.alfa(uu)=alfa;
    set.lambdaIt(uu)=set.lambda;
    set.kappaIt(uu)=set.kappa;  
    minro=min(1./base.m);
    maxro=max(1./base.m);

    

    zero=sparse(base.totP,base.totC); 
    katN=nnz(imag(base.d))+nnz(real(base.d));


    %% each iteration trade-off cooling

    % % if(uu>1 && set.rdmd(uu)>-5)
    % % kat=0.1
    % % else
    % kat=1;
    % % end
    % 
    % set.lambdaIt(uu)=set.lambda*kat;
    % set.lambda=set.lambda*kat;
% 
%     base.CC=[set.lambda*base.C1'*base.C1 zero; zero' set.kappa*base.C2]; 
%     ek1=speye(base.totP)*set.preDiagAdd;
%     ek2=speye(base.totC)*set.preDiagAdd;
%     base.M=[set.lambda*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)]; 
% %     ind=find(diag(base.M)~=0);
% %     ekle=speye(size(base.M,1));
% %     ekle(ind,ind)=0;
% %     base.M=base.M+ekle;
%     base.M=base.ME+base.M;

% ,uu,set.MISFIT(uu)/(length(base.d)*2),...

    
fprintf(['\nINVERSION STEP = %d\n',...
         'misfit = %.2f\n',...
         'rms(all) = %.2f\n',...
         'rdmd(%%) = %.2f\n',...
         'norm(dp) = %.2f\n',...
         'Model Roughness = %.2f\n',...
         'Distortion Strength = %.2f\n',...
         'alfa (step length) = %.5f\n',...
         'min(rho)=%.2e max(rho)=%.2e\n\n']...
         ,uu,set.MISFIT(uu),...         
         set.RMSall(uu),...
         set.rdmd(uu),set.ndp(uu),set.modelr(uu),set.distr(uu),set.alfa(uu),minro,maxro);



lambdaok=0;

if((set.rdmd(uu))>set.limitkappa && uu>1 && (set.sw(5)=='1' || set.sw(9)=='1'))
% if(set.rdmd(uu)>set.limitkappa && uu>1)
    lambdaok=1;   
    zero=sparse(base.totP,base.totC);

    set.kappaIt(uu)=set.kappa/set.kappabol;
    set.kappa=set.kappa/set.kappabol;
    
%     set.lambdaIt(uu)=set.lambda/set.lambdabol;
%     set.lambda=set.lambda/set.lambdabol;

    base.CC=[set.lambda*base.C1'*base.C1 zero; zero' set.kappa*base.C2]; 
    base.CS=[sqrt(set.lambda)*base.C1 zero; zero' sqrt(set.kappa)*base.C2]; 
    
    ek1=speye(base.totP)*set.preDiagAdd;
    ek2=speye(base.totC)*set.preDiagAdd;


    base.M=[set.lambda*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)];

     % base.M=[set.lambdaIt(1)*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)];
   
%     ind=find(diag(base.M)~=0);
%     ekle=speye(size(base.M,1));
%     ekle(ind,ind)=0;
%     base.M=base.M+ekle;

    % base.M=base.ME+base.M;
    % U=ichol(base.M);
    % base.ML=transpose(U);
    % base.MU=U;    

    % [ML,MU]=ilu(base.M);
    % base.ML=ML;
    % base.MU=MU;

    ops.shape='lower';
    base.M=base.ME+base.M;    
    L=ichol(base.M,ops);
    base.MU=transpose(L);
    base.ML=L;      


    
%     if(strcmp(set.incompleteNE,'ilu'))
%     [ML,MU]=ilu(base.M);
%     base.ML=ML;
%     base.MU=MU;
%     elseif(strcmp(set.incompleteNE,'ichol'))
%     U=ichol(base.M);
%     base.ML=transpose(U);
%     base.MU=U;
%     elseif(strcmp(set.incompleteNE,'chol'))
%     U=chol(base.M);
%     base.ML=transpose(U);
%     base.MU=U;
%     else
%     error('Hatalı ichol ya da ilu\n');
%     end
%     
    % set.rs=0;
    set.kappaC=set.kappaC+1;  
  
    fprintf('Kappa is halved\n');
    bas=1;
end

if(set.ndp(uu)<set.limitndp || ne>set.GNlimit)
% if(set.ndp(uu)<set.limitndp || (set.rdmd(uu))>set.limitkappa)
    
    zero=sparse(base.totP,base.totC);


    % if(set.kappaPlambda==1 && lambdaok==0)
    % set.kappaIt(uu)=set.kappa/set.kappabol;
    % set.kappa=set.kappa/set.kappabol;
    % end



    set.lambdaIt(uu)=set.lambda/set.lambdabol;
    set.lambda=set.lambda/set.lambdabol;
    
    base.CC=[set.lambda*base.C1'*base.C1 zero; zero' set.kappa*base.C2];
    base.CS=[sqrt(set.lambda)*base.C1 zero; zero' sqrt(set.kappa)*base.C2]; 
    
    ek1=speye(base.totP)*set.preDiagAdd;
    ek2=speye(base.totC)*set.preDiagAdd;


    base.M=[set.lambda*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)]; 

    % base.M=[set.lambdaIt(1)*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)]; 


%     ind=find(diag(base.M)~=0);
%     ekle=speye(size(base.M,1));
%     ekle(ind,ind)=0;
%     base.M=base.M+ekle;

    % base.M=base.ME+base.M;
    % U=ichol(base.M);
    % base.ML=transpose(U);
    % base.MU=U;  

    ops.shape='lower';
    base.M=base.ME+base.M;    
    L=ichol(base.M,ops);
    base.MU=transpose(L);
    base.ML=L;  



    % [ML,MU]=ilu(base.M);
    % base.ML=ML;
    % base.MU=MU;    

    
%     if(strcmp(set.incompleteNE,'ilu'))
%     [ML,MU]=ilu(base.M);
%     base.ML=ML;
%     base.MU=MU;
%     elseif(strcmp(set.incompleteNE,'ichol'))
%     U=ichol(base.M);
%     base.ML=transpose(U);
%     base.MU=U;
%     else
%     error('Invalid ichol or ilu\n');    
%     end  
    
    % set.rs=0;    
    
    set.ndpC=set.ndpC+1;
    bas=1;    
    fprintf('Norm(dp) is under the specified limit or GN precision is too low, relres=%e \n',ne);
    
    if(set.rdmd(uu)>set.limitkappa && uu>1 && set.sw(5)=='0' && set.sw(9)=='0')    
    set.rs=0;
    set.kappaC=set.kappaC+1;  
    fprintf('kappa counter is increased\n');
    end    
% elseif(set.sw(5)=='1' && set.rdmd(uu)>set.limitkappa && uu>1)
elseif(set.sw(5)=='0' && set.rdmd(uu)>set.limitkappa && uu>1)

    % zero=sparse(base.totP,base.totC);
    % 
    % set.lambdaIt(uu)=set.lambda/set.lambdabol;
    % set.lambda=set.lambda/set.lambdabol;
    % 
    % base.CC=[set.lambda*base.C1'*base.C1 zero; zero' set.kappa*base.C2];
    % base.CS=[sqrt(set.lambda)*base.C1 zero; zero' sqrt(set.kappa)*base.C2]; 
    % 
    % ek1=speye(base.totP)*set.preDiagAdd;
    % ek2=speye(base.totC)*set.preDiagAdd;
    % 
    % 
    % base.M=[set.lambda*(base.C1'*base.C1+ek1) zero; zero' set.kappa*(base.C2+ek2)]; 
    % 
    % base.M=base.ME+base.M;
    % 
    % U=ichol(base.M);
    % base.ML=transpose(U);
    % base.MU=U;    
    % 
    % 
    % set.ndpC=set.ndpC+1;
    % bas=1;    
    % fprintf('No Distortion and RDMD is under limit and lambda is reduced\n');

    fprintf('No Distortion and RDMD is under limit NO action is taken\n');


end


if(bas==1)
    fprintf('Halvings -> kappa=%d\\%d and norm(dp)=%d\\%d \n\n',set.kappaC,set.stopkappa,set.ndpC,set.stopndp);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(set.kappaC>=set.stopkappa)
    set.stop=1;
    fprintf('kappa is halved and the limit is reached \n');    
end

if(set.ndpC>=set.stopndp)
    set.stop=1;
    fprintf('Norm(dp) is under the limit and the limit is reached \n');    
end
    
if(set.RMSall(uu)<=set.stoprms)
    set.stop=1;
    fprintf('RMS is under the limit\n');    
end

if(set.MISFIT(uu)<=set.stopmisfit)
    set.stop=1;
    fprintf('Misfit is under the limit\n');    
end

if(uu>=set.maxitInv)
    set.stop=1;
    fprintf('Maximum number of iterations is reached\n');    
end

n=set.stopposrdmd-1;
st=max(1,uu-n);
en=uu;
ch=sum(set.rdmd(st:en)>0);
if(ch==set.stopposrdmd)
    set.stop=1;
    fprintf('Positive rdmd in the last %d iterations\n',set.stopposrdmd);  
end



end

