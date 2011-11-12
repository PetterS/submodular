function plot_batchrun_new(filename)
    data = load(filename);
    disp(filename)
    

	fprintf('Number of experiments: %d\n',size(data,1));
	
    hocr       = data(:,3+1);
    hocr_itr   = data(:,4+1);
    optimal    = data(:,5+1);
    heuristic  = data(:,6+1);
    fixetal    = data(:,15+1);
    fixetal_itr= data(:,16+1);

    hocrbound      = data(:,7+1);
    optimalbound   = data(:,9+1);
    heuristicbound = data(:,10+1);
    fixetalbound   = data(:,17+1);
    
    disp('HOCR');
    print_bound(optimalbound, hocrbound);
    disp('Fix et al.');
    print_bound(optimalbound, fixetalbound);
    disp('Heuristic');
    print_bound(optimalbound, heuristicbound);
    
    
    edges = linspace(0,data(1,1),20);

    clf;
    hold on
    
    h1 = plot_hist(hocr,      edges,[1 0 0]);
    h2 = plot_hist(fixetal,   edges,[0 1 0]);
    h3 = plot_hist(heuristic, edges,[0 0 1]);
    h4 = plot_hist(optimal,   edges,[0 0 0]);
    
    xlim([0,data(1,1)]);
    
    legend([h1 h2 h3 h4], {'HOCR','Fix et al.','Heuristic','Optimal'});
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

function h = plot_hist(lab,edges,color)
    N1 = histc(lab,edges);
    h = bar(edges,N1,'histc');
    set(h,'FaceColor','none');
    h = stairs(edges,N1,'LineWidth',3,'Color',color);
end