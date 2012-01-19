function plot_batchrun_new(filename,nedges)
    if nargin<2
        nedges=20;
    end

    data = load(filename);
    disp(filename)
    

	fprintf('Number of experiments: %d\n',size(data,1));
	
    hocr       = data(:,3+1);
    hocr_itr   = data(:,4+1);
    optimal    = data(:,5+1);
    heuristic  = data(:,6+1);
    fixetal    = data(:,15+1);
    fixetal_itr= data(:,16+1);
	generators = data(:,21+1);

    hocrbound      = data(:,7+1);
    optimalbound   = data(:,9+1);
    heuristicbound = data(:,10+1);
    fixetalbound   = data(:,17+1);
	generatorsbound= data(:,22+1);
    
	disp('Generators');
    print_bound(optimalbound, generatorsbound);
    disp('HOCR');
    print_bound(optimalbound, hocrbound);
    disp('Fix et al.');
    print_bound(optimalbound, fixetalbound);
    disp('Heuristic');
    print_bound(optimalbound, heuristicbound);
    
    
    edges = linspace(0,data(1,1),nedges);

    clf;
    hold on
    
%     c = colorspiral;
    h1 = plot_hist(hocr,      edges,[ 0.9526    0.5773    0.8113],'-');
    h2 = plot_hist(fixetal,   edges,[ 0.1722    0.7860    0.7948],'-');
    h3 = plot_hist(heuristic, edges,[ 0.6857    0.4521    0.0270],'-');
    h4 = plot_hist(optimal,   edges,[0 0 0],'-');
	h5 = plot_hist(generators,   edges,[0 0 0],':');
    
%     xlim([0,data(1,1)]);
    maxxval = max([hocr(:);hocr_itr(:);optimal(:);heuristic(:);fixetal(:);fixetal_itr(:);generators(:)]);
	xlim([0,1.1*maxxval]);
	
	yl = ylim;
	yl(2) = yl(2)+1;
	ylim(yl);
	
	goodnumberofbins = 50 * 370 / maxxval 
	
%     leg = legend([h1 h2 h3 h4 h5], {'HOCR','Fix et al.','GRD-heuristic','GRD','GRD-Gen.'});
	leg = legend([h5 h4 h3 h2 h1], {'GRD-gen','GRD','GRD-heuristic','Fix et al.','HOCR'});
	set(leg, 'box', 'off');
    xlabel('Number of persistencies');
    ylabel('Frequency');
   
end

function print_bound(opt, bnd)
    relbound = (opt - bnd)./abs(opt);
    minrelbound    = min(relbound);
    medianrelbound = median(relbound);
    maxrelbound    = max(relbound);
    fprintf('min/med/max : %f  %f  %f\n',minrelbound,medianrelbound,maxrelbound);
end

function h = plot_hist(lab,edges,color,type)
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