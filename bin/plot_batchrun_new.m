function plot_batchrun_new(filename,nedges,use_packing)
    if nargin<2
        nedges=20;
	end
	if nargin<3
		use_packing = false;
	end

    data = load(filename);
	if size(data,2) < 24
		data(end,24)=0;
	end
		
    disp(filename)
    

	fprintf('Number of experiments: %d\n',size(data,1));
	
    hocr       = data(:,3+1);
    hocr_itr   = data(:,4+1);
    optimal    = data(:,5+1);
    heuristic  = data(:,6+1);
    fixetal    = data(:,15+1);
    fixetal_itr= data(:,16+1);
	generators = data(:,21+1);
	if use_packing
		packing    = data(:,24+1);
	end

    hocrbound      = data(:,7+1);
    optimalbound   = data(:,9+1);
    heuristicbound = data(:,10+1);
    fixetalbound   = data(:,17+1);
	generatorsbound= data(:,22+1);
	if use_packing
		packingbound = data(:,25+1);
	end
	
	hocrtime        = data(:,11+1);
	optimaltime     = data(:,13+1);
	heuristictime   = data(:,14+1);
	fixetaltime     = data(:,19+1);
	generatorstime  = data(:,23+1);
	if use_packing
		packingtime = data(:,26+1);
	end
	
	%Exclude experiments
% 	heuristic(:) = -1;
    
	disp('Generators');
    print_bound(optimalbound, generatorsbound);
	print_time(generatorstime);
	disp('Optimal');
	print_time(optimaltime);
	disp('Heuristic');
    print_bound(optimalbound, heuristicbound);
	print_time(heuristictime);
	disp('Fix et al.');
    print_bound(optimalbound, fixetalbound);
	print_time(fixetaltime);
    disp('HOCR');
    print_bound(optimalbound, hocrbound);
	print_time(hocrtime);
	if use_packing
		disp('Packing');
		print_bound(optimalbound, packingbound);
		print_time(packingtime);
	end
    
    
    
    edges = linspace(0,data(1,1),nedges);

    clf;
    hold on
    
%     c = colorspiral;
	if ~use_packing
		h1 = plot_hist(hocr,      edges,[ 0.9526    0.5773    0.8113],'-');
		h2 = plot_hist(fixetal,   edges,[ 0.1722    0.7860    0.7948],'-');
		h3 = plot_hist(heuristic, edges,[ 0.6857    0.4521    0.0270],'-');
		h4 = plot_hist(optimal,   edges,[0 0 0],'-');
		h5 = plot_hist(generators,   edges,[0 0 0],':');
		h = [h5 h4 h3 h2 h1];
		leg_str = {'GRD-gen','GRD','GRD-heur.','Fix et al.','HOCR'};
	else
		h1 = plot_hist(hocr,      edges,[ 0.9526    0.5773    0.8113],'-');
		h5 = plot_hist(generators,   edges,[0 0 0],':');
		h6 = plot_hist(packing, edges, [ 0.3126    0.6160    0.9537], '-');
		h = [h5 h6 h1];
		leg_str = {'GRD-gen','Vertex packing','HOCR'};
	end
	
    
    maxxval = max([hocr(:);hocr_itr(:);optimal(:);heuristic(:);fixetal(:);fixetal_itr(:);generators(:)]);
	
	xlim([0,data(1,1)]);
% 	xlim([0,1.05*maxxval]);
	
	yl = ylim;
	yl(2) = yl(2)+1;
	ylim(yl);
	
	goodnumberofbins = 50 * 370 / maxxval 
	
	leg_str = leg_str(h>0);
	h = h(h>0);
	leg = legend(h, leg_str);
	set(leg, 'box', 'off');
    xlabel('Number of persistencies');
    ylabel('Frequency');
   
end

function print_bound(opt, bnd)
	if all(opt>=0) 
		return
	end
    relbound = (opt - bnd)./abs(opt);
    minrelbound    = min(relbound);
    medianrelbound = median(relbound);
    maxrelbound    = max(relbound);
    fprintf('min/med/max : %f  %f  %f\n',minrelbound,medianrelbound,maxrelbound);
end

function print_time(time)
    fprintf('min/med/max : %f  %f  %f  (time)\n',min(time),median(time),max(time));
end

function h = plot_hist(lab,edges,color,type)
	if all(lab<=0) 
		h = 0;
		return
	end
	
    N1 = histc(lab,edges);
    h = bar(edges,N1,'histc');
    set(h,'FaceColor',color);
	set(h,'FaceAlpha',0.3);
    N2 = N1;
    for i = 2:length(N1)-1
        if N1(i-1)==0 && N1(i)==0 && N1(i+1)==0
            N2(i)=nan;
        end
    end
    h = stairs(edges,N2,type,'LineWidth',3,'Color',color);
end