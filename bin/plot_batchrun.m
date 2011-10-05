function plot_batchrun(filename)
    data = load(filename);
    disp(filename)
    figure

    qpbo = data(:,3+1);
    lp  = data(:,5+1);

    qpbobound = data(:,7+1);
    lpbound = data(:,9+1);
    
    relbound = (lpbound - qpbobound)./abs(lpbound);
    minrelbound = min(relbound)
    medianrelbound = median(relbound)
    maxrelbound = max(relbound)


    edges = linspace(0,data(1,1),20);
    NQ = histc(qpbo,edges);
    NL = histc(lp,edges);


    % mx = max(max(qpbo,lp));
    % xmax = 10*ceil(mx/10);
    xmax = data(1,1);

    edges = linspace(0,xmax,20);
    NQ = histc(qpbo,edges);
    NL = histc(lp,edges);

    ymax = max(max(NQ,NL));


    hold on
    hQ = bar(edges,NQ,'histc');
    hL = bar(edges,NL,'histc');
    NN = NQ;
    NN(NQ>NL)=0;
    hM = bar(edges,NN,'histc');
    set(hQ,'FaceColor',[0 0 0]);
    set(hQ,'EdgeColor',[1 1 1]);
    set(hL,'FaceColor',[0.7 0.7 0.7]);
    set(hL,'EdgeColor',[0 0 0]);
    set(hM,'FaceColor',[0 0 0]);
    set(hM,'EdgeColor',[1 1 1]);
    N0=0*NN;
    h0 = bar(edges,N0,'histc');
    set(h0,'FaceColor',[0 0 0]);
    set(h0,'EdgeColor',[0 0 0]);
    xlabel('Number of persistencies');
    ylabel('Frequency');
    xlim([0 xmax]);

    ylim([0 ymax])

    legend({'HOCR','Optimal'})
end