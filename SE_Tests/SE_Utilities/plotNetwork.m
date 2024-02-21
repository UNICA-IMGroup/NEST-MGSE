function [h, p] = plotNetwork(branch_from, branch_to, branch_num, monitoredNodes, monitoredBranches)

h = figure; 

[p, G] = plotGraph(branch_from, branch_to, branch_num);

if ~exist('monitoredNodes') 
    monitoredNodes = [];
end
if ~exist('monitoredBranches')
    monitoredBranches = [];
end
    
%highlight(p, Vmag);
misIndex = ismember(branch_num, monitoredBranches);

title_str = strcat('Nodes: ', num2str(max(branch_to)), '; branches: ', num2str(length(branch_to)));


nVmag = length(monitoredNodes);
nImag_br = size(monitoredBranches, 1);
nMis = nVmag + nImag_br; 

if nMis
    subtitle_str = strcat('measurements points: ', num2str(nMis));

    if nVmag 
        subtitle_str = strcat(subtitle_str, '; nodes: ', num2str(nVmag));
    end
    if nImag_br 
        subtitle_str = strcat(subtitle_str, '; branches: ', num2str(nImag_br));
    end
    
    addMeasurements(p, monitoredNodes, monitoredBranches(:, 1), monitoredBranches(:, 2));
    title({title_str, subtitle_str});
else
    title(title_str);
end
end

function [p, G] = plotGraph(branch_from, branch_to, branch_num, layoutType)

    if nargin < 4
       layoutType = 'layered';
    end
    %G = digraph(branch_from, branch_to);
    G = graph(branch_from, branch_to);
    Labels = cellstr(num2str(branch_num(:)));
    p = plot(G, 'Layout', layoutType);
    
    labeledge(p, branch_from, branch_to, Labels);
    
    p.NodeColor = 'red';
    p.NodeFontSize = 9;
    p.NodeLabelColor = 'b';
    
    p.LineWidth = 1;
    p.EdgeColor = 'k';
    p.EdgeLabelColor = 'm';
    p.EdgeFontSize = 8;
end

function addMeasurements(p, nodi, from, to, cN, cB)

    if nargin < 5
        cN = 'g';
        cB = 'g';
    end
    highlight(p, nodi, 'NodeColor',cN, 'MarkerSize', 7);
    highlight(p, from, to, 'EdgeColor', cB, 'LineWidth', 2);

end