% GRD Matlab wrapper
%
% Petter Strandmark 2012
% petter@maths.lth.se
function varargout = grd(UE,PI,PE,TI,TE,QI,QE,varargin)
	type = class(UE);
	if nargin < 3
		PI = zeros(2,0,'uint32');
		PE = zeros(4,0,type);
	end
	if nargin < 5
		TI = zeros(3,0,'uint32');
		TE = zeros(8,0,type);
	end
	if nargin < 7
		QI = zeros(4,0,'uint32');
		QE = zeros(16,0,type);
	end
	[varargout{1:nargout}] = grd_mex(UE,cast(PI,'uint32'),cast(PE,type), ...
                                     cast(TI,'uint32') ,cast(TE,type), ....
                                     cast(QI,'uint32'), cast(QE,type), ...
                                     varargin{:});
end
