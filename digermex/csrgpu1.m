function [row,col,val,m,n] = csrgpu1(WM)

[val,row,col]=sparse2csr((WM),0);

m=int32((size(WM,1)));
n=int32((size(WM,2)));


val=gpuArray(val);
row=gpuArray(row);
col=gpuArray(col);
end