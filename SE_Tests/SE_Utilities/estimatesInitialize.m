
function estimates = estimatesInitialize(num_nodes, num_branches, num_phases, num_methods, MC)

estimates.V = zeros(num_nodes * num_phases, num_methods, MC);
estimates.I = zeros(num_branches * num_phases, num_methods, MC);
estimates.Time = zeros(num_methods, MC);
%estimates.Iter = zeros(length(SE_type_Group), MC);