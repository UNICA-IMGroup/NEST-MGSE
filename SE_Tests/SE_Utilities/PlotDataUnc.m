
function PlotDataUnc(xData, yData, E, title_txt, xLabel_txt, yLabel_txt, xTicks_data)
    xUnc = [xData, fliplr(xData)];
    yE = E;
    yUnc = [yData + yE, fliplr(yData - yE)];
 
    hold on
    plot(xData, yData, '*b');
    %hPatch = patch(xUnc, yUnc, [0.6 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    fill(xUnc, yUnc, [0.6 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
   
    title(title_txt);
    xticks(xTicks_data);
    xticklabels(xTicks_data);  
    xlabel(xLabel_txt);
    ylabel(yLabel_txt)
    legend ({ 'estimate', 'uncertainty'}, 'Location', 'Best');
    hold off
end