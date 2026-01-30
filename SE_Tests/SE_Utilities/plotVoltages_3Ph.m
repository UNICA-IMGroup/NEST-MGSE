function h = plotVoltages_3Ph(V_true, V_est, topology, SE_type_Group, MC, mciter_to_plot)


mciter_to_plot = min(mciter_to_plot, MC);
nMethods = length(SE_type_Group);

if nMethods > 2
    disp([mfilename, ': alert, this plot routing is customized to work with two methods']);
    return;
end

Vmodule_est = abs(V_est);
Vmodule_true = abs(V_true);
Vangle_est = angle(V_est);
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

Vmodule_est_method1 =  squeeze(Vmodule_est(estimatedNodes, 1, :));

Vmodule_true = Vmodule_true(estimatedNodes, :);

Vangle_est_method1 =  squeeze(Vangle_est(estimatedNodes, 1, :));

Vangle_true = Vangle_true(estimatedNodes, :);

if nMethods > 1
    Vmodule_est_method2 =  squeeze(Vmodule_est(estimatedNodes, 2, :));
    Vangle_est_method2 =  squeeze(Vangle_est(estimatedNodes, 2, :));
end

%% Plot V module estimation, specific iteration
h(1) = figure;
subplot(3, 1, 1)
hold on
plot(nodesA, Vmodule_est_method1(1 : nA, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesA, Vmodule_est_method2(1 : nA, mciter_to_plot), 'sg');
end
plot(nodesA, Vmodule_true(1 : nA, mciter_to_plot), 'Or');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase A'));
ylabel('[pu]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vmodule_est_method1(nA + 1 : nA + nB, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesB, Vmodule_est_method2(nA + 1 : nA + nB, mciter_to_plot), 'sg');
end
plot(nodesB, Vmodule_true(nA + 1 : nA + nB, mciter_to_plot), 'Or');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase B'));
ylabel('[pu]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vmodule_est_method1(nA + nB + 1 : end, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesC, Vmodule_est_method2(nA + nB + 1 : end, mciter_to_plot), 'sg');
end
plot(nodesC, Vmodule_true(nA + nB + 1 : end, mciter_to_plot), 'Or');

xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase C'));
ylabel('[pu]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off
sgtitle(strcat('V module estimation ', ' MC iteration ', num2str(mciter_to_plot)));

%% Plot V angle estimation, specific iteration
h(2) = figure;
subplot(3, 1, 1)
hold on
plot(nodesA, Vangle_est_method1(1 : nA, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesA, Vangle_est_method2(1 : nA, mciter_to_plot), 'sg');
end
plot(nodesA, Vangle_true(1 : nA, mciter_to_plot), 'Or');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase A'));
ylabel('[rad]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vangle_est_method1(nA + 1 : nA + nB, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesB, Vangle_est_method2(nA + 1 : nA + nB, mciter_to_plot), 'sg');
end
plot(nodesB, Vangle_true(nA + 1 : nA + nB, mciter_to_plot), 'Or');
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase B'));
ylabel('[rad]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vangle_est_method1(nA + nB + 1 : end, mciter_to_plot), '*b');
if nMethods > 1
    plot(nodesC, Vangle_est_method2(nA + nB + 1 : end, mciter_to_plot), 'sg');
end
plot(nodesC, Vangle_true(nA + nB + 1 : end, mciter_to_plot), 'Or');

xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase C'));
ylabel('[rad]');
xlabel('nodes');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off
sgtitle(strcat('V angle estimation ', ' MC iteration ', num2str(mciter_to_plot)));
%% plot V module estimation results

Vm_rmse_1 = 100 * rmse_rel(Vmodule_true,  Vmodule_est_method1);
if nMethods > 1
    Vm_rmse_2 = 100 * rmse_rel(Vmodule_true,  Vmodule_est_method2);
end

h(3) = figure;

subplot(3, 1, 1)
hold on
plot(nodesA, Vm_rmse_1(1 : nA), '*b');
if nMethods > 1
    plot(nodesA, Vm_rmse_2(1 : nA), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase A'));
ylabel('RMSE [%]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vm_rmse_1(nA + 1 : nA + nB), '*b');
if nMethods > 1
    plot(nodesB, Vm_rmse_2(nA + 1 : nA + nB), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase B'));
ylabel('RMSE [%]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vm_rmse_1(nA + nB + 1 : end), '*b');
if nMethods > 1
    plot(nodesC, Vm_rmse_2(nA + nB + 1 : end), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
title(strcat('Phase C'));
ylabel('RMSE [%]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off
sgtitle('V module error');

%% plot V angle estimation results
Vangle_rmse_1 = 100 * rmse(Vangle_true,  Vangle_est_method1);
if nMethods > 1
    Vangle_rmse_2 = 100 * rmse(Vangle_true,  Vangle_est_method2);
end

h(4) = figure;
subplot(3, 1, 1)
hold on
plot(nodesA, Vangle_rmse_1(1 : nA), '*b');
if nMethods > 1
    plot(nodesA, Vangle_rmse_2(1 : nA), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase A'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(nodesB, Vangle_rmse_1(nA + 1 : nA + nB), '*b');
if nMethods > 1
    plot(nodesB, Vangle_rmse_1(nA + 1 : nA + nB), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase B'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(nodesC, Vangle_rmse_1(nA + nB + 1 : end), '*b');
if nMethods > 1
    plot(nodesC, Vangle_rmse_1(nA + nB + 1 : end), 'sg');
end
xticks(1 : nNodes);
xticklabels(1 : nNodes);
% xticks(1:length(estimatedNodes));
% xticklabels(estimatedNodes);
title(strcat('Phase C'));
ylabel('RMSE [crad]');
xlabel('nodes');
legend (SE_type_Group, 'Location', 'Best');
hold off
sgtitle('V angle error');

end

function value = rmse(true_value, estimated_value)
    value = sqrt(mean((estimated_value - true_value).^2, 2));
end

function value = rmse_rel(true_value, estimated_value)

value = sqrt(mean(((estimated_value - true_value)./true_value).^2, 2));
end