function measured_values = AddMeasurementErrors_SE(reference_values, Measurements_param)

    measured_values.V = NaN(size(reference_values.V));
    measured_values.I = NaN(size(reference_values.I));
    measured_values.S_br = NaN(size(reference_values.S_br));
    measured_values.S_inj = NaN(size(reference_values.S_inj));


    [std_rel_Vmag_nod, std_rel_Imag_br, ~, ~, ~, ~, ...
        std_rel_PQ_br, std_rel_PQ_inj, ~] = extract_meas_uncertainty(Measurements_param.meas_unc);

    if ~isempty(Measurements_param.meas_indices.Vmag_nod_idx)
        measured_values.V(Measurements_param.meas_indices.Vmag_nod_idx, :) ...
            = abs(reference_values.V(Measurements_param.meas_indices.Vmag_nod_idx, :)) .* (1 + std_rel_Vmag_nod .* randn(size(Measurements_param.meas_indices.Vmag_nod_idx)));
    end
    if ~isempty(Measurements_param.meas_indices.Imag_br_idx)
        measured_values.I(Measurements_param.meas_indices.Imag_br_idx, :) ...
            = abs(reference_values.I(Measurements_param.meas_indices.Imag_br_idx, :)) .* (1 + std_rel_Imag_br .* randn(size(Measurements_param.meas_indices.Imag_br_idx)));
    end
    if ~isempty(Measurements_param.meas_indices.PQ_br_idx)
        measured_values.S_br(Measurements_param.meas_indices.PQ_br_idx, :) ...
            = reference_values.S_br(Measurements_param.meas_indices.PQ_br_idx, :) .* (1 + std_rel_PQ_br .* randn(size(Measurements_param.meas_indices.PQ_br_idx)));
    end
    if ~isempty(Measurements_param.meas_indices.PQ_inj_idx)
        measured_values.S_inj(Measurements_param.meas_indices.PQ_inj_idx, :) ...
            = reference_values.S_inj(Measurements_param.meas_indices.PQ_inj_idx, :) .* (1 + std_rel_PQ_inj .* randn(size(Measurements_param.meas_indices.PQ_inj_idx)));
    end
 end

