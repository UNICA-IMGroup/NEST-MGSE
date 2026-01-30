function h = plotCurrents_3Ph_Advanced(I_true, I_est, topology, MC, mciter_to_plot, bStatistics)
if ~exist('bStatistics', 'var')
    bStatistics = false;
end

mciter_to_plot = min(mciter_to_plot, MC);

Imodule_true = abs(I_true);
Iangle_true = angle(I_true);

nBranches = length(Imodule_true(:, 1)) / 3;
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
Imodule_true = Imodule_true(estimatedBranches, :);

Iangle_est =  I_est.Angle(estimatedBranches, :);
Iangle_est_Unc =  I_est.Unc.Angle(estimatedBranches, :);
Iangle_true = Iangle_true(estimatedBranches, :);


%% Plot I module estimation, specific iteration
h(1) = figure;
yLabel_txt = '[pu]';
xLabel_txt = 'branches';
xTicks_data = 1 : nBranches;

subplot(3, 1, 1)
Data_Indeces = 1 : nA;
xData = reshape(branchesA, 1, []);
yData = reshape(Imodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Imodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Module Phase A');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 2)
Data_Indeces = nA + 1 : nA + nB;
xData = reshape(branchesB, 1, []);
yData = reshape(Imodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Imodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Module Phase B');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);


subplot(3, 1, 3)
Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(branchesC, 1, []);
yData = reshape(Imodule_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Imodule_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Imodule_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Module Phase C');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

%% Plot I angle estimation, specific iteration
yLabel_txt = '[rad]';
h(2) = figure;
subplot(3, 1, 1)
Data_Indeces = 1 : nA;
xData = reshape(branchesA, 1, []);
yData = reshape(Iangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Iangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Angle Phase A');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 2)
Data_Indeces = nA + 1 : nA + nB;
xData = reshape(branchesB, 1, []);
yData = reshape(Iangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Iangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Angle Phase B');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

subplot(3, 1, 3)
Data_Indeces = nA + nB + 1 : nA + nB + nC;
xData = reshape(branchesC, 1, []);
yData = reshape(Iangle_est(Data_Indeces, mciter_to_plot), 1, []);
yE    = reshape(Iangle_est_Unc(Data_Indeces, mciter_to_plot), 1, []);
yTrue = Iangle_true(Data_Indeces, mciter_to_plot);
title_txt = strcat('Current Angle Phase C');
PlotDataUncReference(xData, yData, yE, yTrue, title_txt, xLabel_txt, yLabel_txt, xTicks_data);

if ~bStatistics
    return;
end
%% plot I module estimation results

Im_rmse_1 = 100 * rmse_rel(Imodule_true,  Imodule_est);


h(3) = figure;

subplot(3, 1, 1)
hold on
plot(branchesA, Im_rmse_1(1 : nA), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase A'));
ylabel('RMSE [%]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Im_rmse_1(nA + 1 : nA + nB), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase B'));
ylabel('RMSE [%]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Im_rmse_1(nA + nB + 1 : end), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase C'));
ylabel('RMSE [%]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off
sgtitle('I module error');

%% plot I angle estimation results
Iangle_rmse_1 = 100 * rmse(Iangle_true,  Iangle_est);

h(4) = figure;
subplot(3, 1, 1)
hold on
plot(branchesA, Iangle_rmse_1(1 : nA), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
% xticks(1:length(estimatedBranches));
% xticklabels(estimatedBranches);
title(strcat('Phase A'));
ylabel('RMSE [crad]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Iangle_rmse_1(nA + 1 : nA + nB), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
% xticks(1:length(estimatedBranches));
% xticklabels(estimatedBranches);
title(strcat('Phase B'));
ylabel('RMSE [crad]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Iangle_rmse_1(nA + nB + 1 : end), '*b');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
% xticks(1:length(estimatedBranches));
% xticklabels(estimatedBranches);
title(strcat('Phase C'));
ylabel('RMSE [crad]');
xlabel('branches');
legend ('estimate', 'Location', 'Best');
hold off
sgtitle('I angle error');

end

function value = rmse(true_value, estimated_value)
    value = sqrt(mean((estimated_value - true_value).^2, 2));
end

function value = rmse_rel(true_value, estimated_value)

value = sqrt(mean(((estimated_value - true_value)./true_value).^2, 2));
end