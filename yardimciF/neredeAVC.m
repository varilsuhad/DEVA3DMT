function [ say ] = neredeAVC(i,j,k,nx,ny,nz,no )

noAx=nx*(ny-1)*(nz-1);
noAy=(nx-1)*(ny)*(nz-1);
noAz=(nx-1)*(ny-1)*(nz);

switch no
    
    case 1
     if (i>nx || i<1 || j>ny || j<2 || k>nz || k<2)
     fprintf('nx=%d ny%d nz=%d i=%d j=%d k%d\n',nx,ny,nz,i,j,k);
     error('Ax index error\n')            
     else
     xyon=(nz-1)*(ny-1);
     yyon=(nz-1);  
     say=(i-1)*xyon+(j-2)*yyon+(k-1);
     end
        
     case 2
        
     if (i>nx || i<2 || j>ny || j<1 || k>nz || k<2)
         fprintf('%d %d %d\n',i,j,k);
     error('Ay index error\n')            
     else
     xyon=(nz-1)*(ny);
     yyon=(nz-1);  
     say=noAx+(i-2)*xyon+(j-1)*yyon+(k-1);
     end       
           
     case 3
        
     if (i>nx || i<2 || j>ny || j<2 || k>nz || k<1)
     error('Az index error\n')            
     else
     xyon=(nz)*(ny-1);
     yyon=(nz);  
     say=noAx+noAy+(i-2)*xyon+(j-2)*yyon+k;
     end 
     
    case 4
     if (i>nx || i<2 || j>ny || j<2 || k>nz || k<2)
     error('fi index error\n')                
     else
     xyon=(nz-1)*(ny-1);
     yyon=(nz-1);
     say=noAx+noAy+noAz+(i-2)*xyon+(j-2)*yyon+(k-1);    
     end
             
     otherwise
    error('no error\n')            


end

end

