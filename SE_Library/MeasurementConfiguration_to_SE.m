function meas = MeasurementConfiguration_to_SE(num_nodes, num_branches, DG_position, placement_config, unc_config)
% The routine returns the meas struct containg the measurement parameters,
% i.e. the placement index of the intruments and their uncertainty
% Input 
% -num_nodes 
% -num_branches
% -DG_position, vector containing the position of DG
% -placement_config, avilable values {'default', '3 nodes placement'}. 
%  It is possible to add new configurations in the file .Measurements/Measurement_placement.m
% -unc_config, available values {'default', 'nominal'}
%  It is possible to add new configurations in the file .Measurements/Measurement_placement.m
%
% Output
% -meas
%  meas.meas_indices
%  meas.meas_unc
if ~exist('placement_config', 'var')
    placement_config = 'default';
end
if ~exist('unc_config', 'var') 
  unc_config = 'default';  
end

meas.meas_indices = Measurement_placement(placement_config);

% The DG is monitored with magnitude voltage measurement and PQ
% injection
if exist('DG_position', 'var') && ~isempty(DG_position)

    DG_position_toUse = unique(DG_position(DG_position <= num_nodes));
    if DG_position(end) > num_nodes
        disp([mfilename, ': DG positions greater than ', num2str(num_nodes), ' not allowed']);
    end
        
    meas.meas_indices.Vmag_nod_idx = unique([meas.meas_indices.Vmag_nod_idx; DG_position_toUse(:)]);
    meas.meas_indices.PQ_inj_idx = unique([meas.meas_indices.PQ_inj_idx; DG_position_toUse(:)]);
end

% indeces saturation
meas.meas_indices.PQ_br_idx = meas.meas_indices.PQ_br_idx(meas.meas_indices.PQ_br_idx <= num_branches);
meas.meas_indices.Vmag_nod_idx = meas.meas_indices.Vmag_nod_idx(meas.meas_indices.Vmag_nod_idx <= num_nodes);

if isfield(meas.meas_indices, 'Imag_br_idx')
    meas.meas_indices.Imag_br_idx = meas.meas_indices.Imag_br_idx(meas.meas_indices.Imag_br_idx <= num_branches);
end

meas.meas_unc = Measurement_uncertainty(unc_config);

