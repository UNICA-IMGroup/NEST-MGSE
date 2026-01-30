function h = plot_MGSE(V_est, I_est, topology, main_title)

%% Voltages Pretaration
nNodes = topology.num_nodes;

bNodes = reshape(topology.node_3ph, [], 3);
nodesA = find(bNodes(:, 1));               
nA = length(nodesA);

nodesB = find(bNodes(:, 2));  
nB = length(nodesB);

nodesC = find(bNodes(:, 3));  
nC = length(nodesC);
 
estimatedNodes = [nodesA; nNodes + nodesB; nNodes * 2 + nodesC];

Vmodule_est = V_est.Module(estimatedNodes, :);
Vmodule_est_Unc = V_est.Unc.Module(estimatedNodes, :);

Vangle_est = V_est.Angle(estimatedNodes, :);
Vangle_est_Unc = V_est.Unc.Angle(estimatedNodes, :);

%% Plot V module estimation, specific iteration
h = figure;

monitors = get(0, 'MonitorPositions'); % Obtain info on all monitors
        
% Select the monitor desiderato (e.g. the furthest to the left)
[~, idxLeft] = min(monitors(:, 1));
monitorLeft = monitors(idxLeft, :);  % [left, bottom, width, height]

% Imposta la figura con dimensioni in pixel su quel monitor
%figure(h);
set(h, 'Units', 'pixels', 'Position', monitorLeft);
sgtitle(main_title);

yLabel_txt = '[pu]';
xLabel_txt = 'nodes';
xTicks_data = 1 : nNodes;

subplot(3, 4, 1)

Data_Indeces = 1 : nA;
xData = reshape(nodesA, 1, []);
yData = reshape(Vmodule_est(Data_Indeces), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Module Phase A');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 5)

Data_Indeces = nA + 1 : nA + nB;
xData = reshape(nodesB, 1, []);
yData = reshape(Vmodule_est(Data_Indeces), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Module Phase B');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);


subplot(3, 4, 9)

Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(nodesC, 1, []);
yData = reshape(Vmodule_est(Data_Indeces), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Module Phase C');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);



%% Plot V angle estimation, specific iteration
yLabel_txt = '[rad]';
%h(2) = figure;

subplot(3, 4, 2)

Data_Indeces = 1 : nA;
xData = reshape(nodesA, 1, []);
yData = reshape(Vangle_est(Data_Indeces), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Angle Phase A');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 6)

Data_Indeces = nA + 1 : nA + nB;
xData = reshape(nodesB, 1, []);
yData = reshape(Vangle_est(Data_Indeces), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Angle Phase B');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 10)

Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(nodesC, 1, []);
yData = reshape(Vangle_est(Data_Indeces), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Voltage Angle Phase C');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

%% Currents Preparation
%%

nBranches = topology.num_branches;
bBranches = reshape(topology.branch_3ph, [], 3);

branchesA = find(bBranches(:, 1));
nA = length(branchesA);
branchesB = find(bBranches(:, 2));
nB = length(branchesB);
branchesC = find(bBranches(:, 3));
nC = length(branchesC);
 
estimatedBranches = [branchesA; nBranches + branchesB; nBranches * 2 + branchesC];

Imodule_est =  I_est.Module(estimatedBranches, :);
Imodule_est_Unc =  I_est.Unc.Module(estimatedBranches, :);

Iangle_est =  I_est.Angle(estimatedBranches, :);
Iangle_est_Unc =  I_est.Unc.Angle(estimatedBranches, :);

%% Plot I module estimation, specific iteration
%h(2) = figure;
yLabel_txt = '[pu]';
xLabel_txt = 'branches';
xTicks_data = 1 : nBranches;

subplot(3, 4, 3)

Data_Indeces = 1 : nA;
xData = reshape(branchesA, 1, []);
yData = reshape(Imodule_est(Data_Indeces), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Module Phase A');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 7)

Data_Indeces = nA + 1 : nA + nB;
xData = reshape(branchesB, 1, []);
yData = reshape(Imodule_est(Data_Indeces), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Module Phase B');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);


subplot(3, 4, 11)

Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(branchesC, 1, []);
yData = reshape(Imodule_est(Data_Indeces), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Module Phase C');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

%% Plot I angle estimation, specific iteration
yLabel_txt = '[rad]';
%h(4) = figure;
subplot(3, 4, 4)

Data_Indeces = 1 : nA;
xData = reshape(branchesA, 1, []);
yData = reshape(Iangle_est(Data_Indeces), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Angle Phase A');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 8)

Data_Indeces = nA + 1 : nA + nB;
xData = reshape(branchesB, 1, []);
yData = reshape(Iangle_est(Data_Indeces), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Angle Phase B');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 4, 12)

Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(branchesC, 1, []);
yData = reshape(Iangle_est(Data_Indeces), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces), 1, []);
title_txt = strcat('Current Angle Phase C');
PlotDataUnc(xData, yData, yE, title_txt, xLabel_txt, yLabel_txt, xTicks_data);


