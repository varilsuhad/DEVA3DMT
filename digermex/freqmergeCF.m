% © 2020–2025 Deniz Varılsüha — Non-commercial research use only. See LICENSE.
% Contact: deniz.varilsuha@itu.edu.tr
function [pase] = freqmergeCF(pase,f,set)


if(set.MatrixStacking==1)

[pase.rowA,pase.colA,pase.valA,pase.rowM,pase.colM,pase.valM]=freqmergemCF(pase.T1r,pase.T1c,pase.T1v,pase.T2r,pase.T2c,pase.T2v,pase.T3r,pase.T3c,pase.T3v...
    ,pase.Q2r,pase.Q2c,pase.Q2v,double(f));
pase.bstack=gpuArray(pase.b(:));
else
[pase.rowA,pase.colA,pase.valA,pase.rowM,pase.colM,pase.valM]=freqmergemSCF(pase.T1r,pase.T1c,pase.T1v,pase.T2r,pase.T2c,pase.T2v,pase.T3r,pase.T3c,pase.T3v...
    ,pase.Q2r,pase.Q2c,pase.Q2v,double(f));

pase.bstack=gpuArray(pase.b(:));
end


end