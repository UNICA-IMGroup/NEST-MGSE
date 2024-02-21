function values = VIS_inizialize_MC(num_nodes, num_branches, num_phases, MC)

values.V = zeros(num_nodes * num_phases, MC);
values.I = zeros(num_branches * num_phases, MC);
values.S_br = zeros(num_branches * num_phases, MC);
values.S_inj = zeros(num_nodes * num_phases, MC);