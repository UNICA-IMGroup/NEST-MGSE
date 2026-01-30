function h = plotVoltages_3Ph_Advanced(V_true, V_est, topology, MC, mciter_to_plot, bStatistics)
if ~exist('bStatistics', 'var')
    bStatistics = false;
end
mciter_to_plot = min(mciter_to_plot, MC);

Vmodule_true = abs(V_true);
Vangle_true = angle(V_true);
nNodes = length(Vmodule_true(:, 1)) / 3;

bNodes = reshape(topology.node_3ph, [], 3);
nodesA = find(bNodes(:, 1));               
nA = length(nodesA);

nodesB = find(bNodes(:, 2));  
nB = length(nodesB);

nodesC = find(bNodes(:, 3));  
nC = length(nodesC);
 
estimatedNodes = [nodesA; nNodes + nodesB; nNodes * 2 + nodesC];

Vmodule_est =  V_est.Module(estimatedNodes, :);
Vmodule_est_Unc = V_est.Unc.Module(estimatedNodes, :);
Vmodule_true = Vmodule_true(estimatedNodes, :);

Vangle_est = V_est.Angle(estimatedNodes, :);
Vangle_est_Unc = V_est.Unc.Angle(estimatedNodes, :);
Vangle_true = Vangle_true(estimatedNodes, :);

%% Plot V module estimation, specific iteration
h(1) = figure;
yLabel_txt = '[pu]';
xLabel_txt = 'nodes';
xTicks_data = 1 : nNodes;

subplot(3, 1, 1)
Data_Indeces = 1 : nA;
xData = reshape(nodesA, 1, []);
yData = reshape(Vmodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vmodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Module Phase A');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 2)
Data_Indeces = nA + 1 : nA + nB;
xData = reshape(nodesB, 1, []);
yData = reshape(Vmodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vmodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Module Phase B');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);


subplot(3, 1, 3)
Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(nodesC, 1, []);
yData = reshape(Vmodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vmodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vmodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Module Phase C');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);



%% Plot V angle estimation, specific iteration
yLabel_txt = '[rad]';
h(2) = figure;

subplot(3, 1, 1)
Data_Indeces = 1 : nA;
xData = reshape(nodesA, 1, []);
yData = reshape(Vangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Angle Phase A');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 2)
Data_Indeces = nA + 1 : nA + nB;
xData = reshape(nodesB, 1, []);
yData = reshape(Vangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Angle Phase B');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 3)
Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(nodesC, 1, []);
yData = reshape(Vangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Vangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Vangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Voltage Angle Phase C');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);
%sgtitle(strcat('V angle estimation ', ' MC iteration ', num2str(mciter_to_plot)));
if ~bStatistics
    return;
end

%% plot V module estimation results

Vm_rmse_1 = 100 * rmse_rel(Vmodule_true,  Vmodule_est);

h(3) = figure;

subplot(3, 1, 1)
hold on
plot(nodesA, Vm_rmse_1(1 : nA), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase A'));
ylabel('RMSE [%]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vm_rmse_1(nA + 1 : nA + nB), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase B'));
ylabel('RMSE [%]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vm_rmse_1(nA + nB + 1 : end), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase C'));
ylabel('RMSE [%]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off
sgtitle('V module error');

%% plot V angle estimation results
Vangle_rmse_1 = 100 * rmse(Vangle_true,  Vangle_est);

h(4) = figure;
subplot(3, 1, 1)
hold on
plot(nodesA, Vangle_rmse_1(1 : nA), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase A'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vangle_rmse_1(nA + 1 : nA + nB), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase B'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vangle_rmse_1(nA + nB + 1 : end), '*b');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase C'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend ('estimate', 'Location', 'Best');
hold off
sgtitle('V angle error');

end

function value = rmse(true_value, estimated_value)
    value = sqrt(mean((estimated_value - true_value).^2, 2));
end

function value = rmse_rel(true_value, estimated_value)

value = sqrt(mean(((estimated_value - true_value)./true_value).^2, 2));
end


