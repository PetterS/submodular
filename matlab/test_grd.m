clear all
clc

if 0
	delete(['grd_mex.' mexext]);
	rehash
end

UE = rand(2,4);

PI = [1 1 1 2 2 3;
      2 3 4 3 4 4];
PE = rand(4,6);

TI = [1 1;
      2 2; 
	  3 4];
TE = 10*rand(8, size(TI,2));

[L1 energy] = grd(UE, PI, PE, TI, TE, [], [], 'method', 'GRD')

PI = uint32(PI);
TI = uint32(TI);
L2 = vgg_qpbo(UE, PI, PE, TI, TE)
