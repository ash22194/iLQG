function plotForwardPass(X, plot_pairs, plot_labels)
%     fig = get(groot,'CurrentFigure');
%     if (isempty(fig))
        figure;
%     end
    
    num_plots = size(plot_pairs,1);
    for n = 1:1:num_plots
        subplot(num_plots, 1, n);
        hold on;
        x = X(plot_pairs(n,1), :);
        y = X(plot_pairs(n,2), :);
        plot(x, y);
        scatter(x(1), y(1), 20, [0,1,0],'filled');
        text(x(1),y(1),'Start');
        scatter(x(end), y(end), 20, [1,0,0],'filled');
        text(x(end),y(end),'End');
        xlabel(plot_labels(n,1));
        ylabel(plot_labels(n,2));
    end
    
end