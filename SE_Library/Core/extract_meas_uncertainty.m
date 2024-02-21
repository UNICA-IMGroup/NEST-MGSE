
function [std_rel_Vmag_nod, std_rel_Imag_br, std_rel_Isync_mag, std_rel_Isync_phase, std_rel_Vsync_mag, std_rel_Vsync_phase, ...
    std_rel_PQ_br, std_rel_PQ_inj, std_rel_PQ_pseudo] = extract_meas_uncertainty(meas_unc, bUniform)

    cov_factor = 3;
    if exist('bUniform', 'var') && bUniform
        cov_factor = sqrt(3);
    end
    
    cov_factor_100 = cov_factor * 100;
    
    std_rel_Vmag_nod = meas_unc.unc_percent_Vmag_nod / cov_factor_100;
    std_rel_Imag_br = meas_unc.unc_percent_Imag_br / cov_factor_100;
    std_rel_Isync_mag = []; %meas_unc.unc_percent_Isync_mag / cov_factor_100;
    std_rel_Isync_phase = []; %meas_unc.unc_crad_Isync_phase / cov_factor_100;
    std_rel_Vsync_mag = []; %meas_unc.unc_percent_Vsync_mag / cov_factor_100;
    std_rel_Vsync_phase = [];%meas_unc.unc_crad_Vsync_phase / cov_factor_100;
    std_rel_PQ_br = meas_unc.unc_percent_PQ_br / cov_factor_100;
    std_rel_PQ_inj = meas_unc.unc_percent_PQ_inj / cov_factor_100;
    std_rel_PQ_pseudo = meas_unc.unc_percent_PQ_pseudo / cov_factor_100;
end