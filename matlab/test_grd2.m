clc

%
%  f(x) = x1*x2*x3*x4
%
UE = zeros(2,4);
QI = [1;2;3;4];
QE = zeros(16,1);
QE(16,1) = 1;

Methods = {'GRD','GRD-heur','Fix','HOCR'};
for method = Methods
	method=method{1};
	
	disp(method);
	[L bound] = grd(UE, [], [], [], [], QI, QE, 'method', method);
	disp(L');
	disp(bound);
end


% 
