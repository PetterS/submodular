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
TE = rand(8, size(TI,2));

[L1 energy] = grd(UE, PI, PE, TI, TE, [], [], 'method', 'GRD')


%
% This is a wrapper for roof duality by VGG, Oxford.
% Available in Oliver Woodford's imrender toolbox
%
L2 = vgg_qpbo(UE, uint32(PI), PE, uint32(TI), TE)
