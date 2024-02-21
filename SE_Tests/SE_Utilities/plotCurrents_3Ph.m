function plotCurrents_3Ph(I_true, I_est, topology, SE_type_Group, MC, mciter_to_plot)


mciter_to_plot = min(mciter_to_plot, MC);
nMethods = length(SE_type_Group);

if nMethods ~= 2
    disp([mfilename, ': alert, this plot routing is customized to work with two methods']);
    return;
end

Imodule_est = abs(I_est);
Imodule_true = abs(I_true);
Iangle_est = angle(I_est);
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

Imodule_est_method1 =  squeeze(Imodule_est(estimatedBranches, 1, :));
Imodule_est_method2 =  squeeze(Imodule_est(estimatedBranches, 2, :));
Imodule_true = Imodule_true(estimatedBranches, :);

Iangle_est_method1 =  squeeze(Iangle_est(estimatedBranches, 1, :));
Iangle_est_method2 =  squeeze(Iangle_est(estimatedBranches, 2, :));
Iangle_true = Iangle_true(estimatedBranches, :);

%% Plot I module estimation, specific iteration
figure;
subplot(3, 1, 1)
hold on
plot(branchesA, Imodule_est_method1(1 : nA, mciter_to_plot), '*b');
plot(branchesA, Imodule_est_method2(1 : nA, mciter_to_plot), 'sg');
plot(branchesA, Imodule_true(1 : nA, mciter_to_plot), 'Or');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase A'));
ylabel('[A]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Imodule_est_method1(nA + 1 : nA + nB, mciter_to_plot), '*b');
plot(branchesB, Imodule_est_method2(nA + 1 : nA + nB, mciter_to_plot), 'sg');
plot(branchesB, Imodule_true(nA + 1 : nA + nB, mciter_to_plot), 'Or');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase B'));
ylabel('[A]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Imodule_est_method1(nA + nB + 1 : end, mciter_to_plot), '*b');
plot(branchesC, Imodule_est_method2(nA + nB + 1 : end, mciter_to_plot), 'sg');
plot(branchesC, Imodule_true(nA + nB + 1 : end, mciter_to_plot), 'Or');

xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase C'));
ylabel('[A]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off
sgtitle(strcat('I module estimation ', ' MC iteration ', num2str(mciter_to_plot)));

%% Plot I angle estimation, specific iteration
figure;
subplot(3, 1, 1)
hold on
plot(branchesA, Iangle_est_method1(1 : nA, mciter_to_plot), '*b');
plot(branchesA, Iangle_est_method2(1 : nA, mciter_to_plot), 'sg');
plot(branchesA, Iangle_true(1 : nA, mciter_to_plot), 'Or');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase A'));
ylabel('[rad]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Iangle_est_method1(nA + 1 : nA + nB, mciter_to_plot), '*b');
plot(branchesB, Iangle_est_method2(nA + 1 : nA + nB, mciter_to_plot), 'sg');
plot(branchesB, Iangle_true(nA + 1 : nA + nB, mciter_to_plot), 'Or');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase B'));
ylabel('[rad]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Iangle_est_method1(nA + nB + 1 : end, mciter_to_plot), '*b');
plot(branchesC, Iangle_est_method2(nA + nB + 1 : end, mciter_to_plot), 'sg');
plot(branchesC, Iangle_true(nA + nB + 1 : end, mciter_to_plot), 'Or');

xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase C'));
ylabel('[rad]');
xlabel('branches');
legend ([SE_type_Group, 'reference'], 'Location', 'Best');
hold off
sgtitle(strcat('I angle estimation ', ' MC iteration ', num2str(mciter_to_plot)));
%% plot I module estimation results

Im_rmse_1 = 100 * rmse_rel(Imodule_true,  Imodule_est_method1);
Im_rmse_2 = 100 * rmse_rel(Imodule_true,  Imodule_est_method2);


figure;

subplot(3, 1, 1)
hold on
plot(branchesA, Im_rmse_1(1 : nA), '*b');
plot(branchesA, Im_rmse_2(1 : nA), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase A'));
ylabel('RMSE [%]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Im_rmse_1(nA + 1 : nA + nB), '*b');
plot(branchesB, Im_rmse_2(nA + 1 : nA + nB), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase B'));
ylabel('RMSE [%]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Im_rmse_1(nA + nB + 1 : end), '*b');
plot(branchesC, Im_rmse_2(nA + nB + 1 : end), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
title(strcat('Phase C'));
ylabel('RMSE [%]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off
sgtitle('I module error');

%% plot I angle estimation results
Iangle_rmse_1 = 100 * rmse(Iangle_true,  Iangle_est_method1);
Iangle_rmse_2 = 100 * rmse(Iangle_true,  Iangle_est_method2);

figure
subplot(3, 1, 1)
hold on
plot(branchesA, Iangle_rmse_1(1 : nA), '*b');
plot(branchesA, Iangle_rmse_2(1 : nA), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
xticks(1:length(estimatedBranches));
xticklabels(estimatedBranches);
title(strcat('Phase A'));
ylabel('RMSE [crad]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 2)
hold on
plot(branchesB, Iangle_rmse_1(nA + 1 : nA + nB), '*b');
plot(branchesB, Iangle_rmse_1(nA + 1 : nA + nB), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
xticks(1:length(estimatedBranches));
xticklabels(estimatedBranches);
title(strcat('Phase B'));
ylabel('RMSE [crad]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off

subplot(3, 1, 3)
hold on
plot(branchesC, Iangle_rmse_1(nA + nB + 1 : end), '*b');
plot(branchesC, Iangle_rmse_1(nA + nB + 1 : end), 'sg');
xticks(1 : nBranches);
xticklabels(1 : nBranches);
xticks(1:length(estimatedBranches));
xticklabels(estimatedBranches);
title(strcat('Phase C'));
ylabel('RMSE [crad]');
xlabel('branches');
legend (SE_type_Group, 'Location', 'Best');
hold off
sgtitle('I angle error');

end

function value = rmse(true_value, estimated_value)
    value = sqrt(mean((estimated_value - true_value).^2, 2));
end

function value = rmse_rel(true_value, estimated_value)

value = sqrt(mean(((estimated_value - true_value)./true_value).^2, 2));
end