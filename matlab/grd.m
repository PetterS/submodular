% GRD Matlab wrapper for generalized roof duality
%
% Using a similar interface as vgg_qpbo by VGG, Oxford
%
%   [L energy] = grd(UE, PI, PE, TI, TE, QI, QE, [options])
%
% IN:
%   UE : 2-by-n matrix of unary terms for M nodes
%        [ E1(0) E2(0) ... EM(0);
%          E1(1) E2(1) ... EM(1)]
%
%   PI : 2-by-P or 3-by-P uint32 matrix. 
%        Each column is [i; j; ind], where i and j are the clique indices
%        and ind is the index to the pairwise energy table. If ind is
%        omitted, the same column in the energy table is used
%   PE : 4-by-p pairwise energy table.
%        Each column is [E00; E01; E10; E11] for a given clique
%
%   TI, TE, QI and QE are used in the same way. 
%		TI : 3-by-T or 4-by-T
%		TE : 8-by-t, each column [E000; E001; .... E111]
%		QI : 4-by-Q or 5-by-Q
%		QE : 16-by-q, each column [E0000; E0001; .... E1111]
%
%   UE can be either double or int32. PE, TE and QE will be cast to the
%   type of UE.
%
%
%  Options are specified as 
%	
%	 [...] = grd(..., 'option', value)
%
%	'method'	:	{'GRD','GRD-gen','GRD-heur','Fix','HOCR'}
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
