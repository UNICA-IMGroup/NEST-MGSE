function [h, nMis]  = PlotNetwork_to_SE(nominal_linedata_t, meas_indices)

    monitoredNodes = unique([meas_indices.PQ_inj_idx; meas_indices.Vmag_nod_idx]);
    monitoredBranches_index = unique([meas_indices.Imag_br_idx; meas_indices.PQ_br_idx]);
    monitoredBranches = [nominal_linedata_t.branch_first(monitoredBranches_index, :), nominal_linedata_t.branch_end(monitoredBranches_index, :)]; 

    [h, p, nMis] = plotNetwork(nominal_linedata_t.branch_first, nominal_linedata_t.branch_end, nominal_linedata_t.branch_cod, monitoredNodes, monitoredBranches);    