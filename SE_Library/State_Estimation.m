function [V, I, num_iter] = State_Estimation(Measurements, pseudomeas, Measurements_param, Network_param, SE_type, bInitialize)
% The routine performes the state estimation, that is the estimation of all
% the N node voltage phasors and all the L branch current phasors 
% input:
% - Measurements: struct of measurement vectors 
%   V voltages, I Currents, S_br power flows, S_inj power injections
% - pseudomeas: struct containing pseudo-measurements vectors
%   P active power, Q rective power
% - Measurements_param describes the measurement system:
%   meas_indices contains the vector of the indices for each type if
%   available measurement, while meas_unc contains the uncertainty of each
%   measurement and pseudo-measurement.
% - Network_param describes the network using the fields: topolgy, Zbr_set
%   (impedeces of the brenches), nominal_linedata (table that collect useful network
%   descriction for the State Estimation), mesh describes the mesh
%   information
% - SE_type has two possible values "BC rect" and "NV rect" for Brench
%   currents and Node Voltages state estimaton
% - bInitialize is an optional boolean parameter to force the
%   initialization of matrices used by the state estimation. 
% output: 
% -V node voltages in real V.Rect, and imaginary part V.Imag. Both of them
%  have N x 3 dimension (assuming zero value when the node is not available
%  for the specific phase).
% -I branch Currents in real I.Rect, and imaginary part I.Imag. Both of them
%  have L x 3 dimension (assuming zero value when the node is not available
%  for the specific phase)


if ~exist('bInitialize', 'var')
    bInitialize = false;
end
persistent YZ_Matrices;
persistent Mesh_Matrices;

if isempty(YZ_Matrices) || bInitialize
    [YZ_Matrices, Mesh_Matrices] = Impedence_Matrices_to_SE(Network_param);
end

num_nodes = Network_param.topology.num_nodes;
num_branches = Network_param.topology.num_branches;
node_3ph = Network_param.topology.node_3ph;
node_3ph_reverse = Network_param.topology.node_3ph_reverse;
nominal_linedata = table2array(Network_param.nominal_linedata);

% it makes the measurements informations available for the state estimation
zdata = Zdata_3Ph_to_SE(nominal_linedata, pseudomeas, Measurements, Measurements_param);

if ~exist('SE_type', 'var')
    SE_type = 'BC_rect';
end

if ~exist('V_abs', 'var')
    V_abs = ones (3 * num_nodes, 1);
end
if ~exist('V_theta', 'var')
    V_theta = [zeros(num_nodes,1); ones(num_nodes, 1) * (2 * pi/3); ones(num_nodes, 1) * (-2 * pi / 3)];
end

switch SE_type
    case 'BC rect'
        [V.Real, V.Imag, I.Real, I.Imag, num_iter] = BC_SE_rect_3Ph(num_nodes, num_branches, nominal_linedata, node_3ph, node_3ph_reverse, V_abs, V_theta, ...
            table2array(zdata), YZ_Matrices.Znod, YZ_Matrices.Znod_state, Mesh_Matrices.Gzcomplex, Mesh_Matrices.Gz, Mesh_Matrices.Gz_state);
        
    case 'NV rect'
        [V.Real, V.Imag, I.Real, I.Imag, num_iter] = NV_SE_rect_3Ph(num_nodes, num_branches, nominal_linedata, node_3ph, node_3ph_reverse, V_abs, V_theta, ...
            table2array(zdata), YZ_Matrices.Y_state, YZ_Matrices.Y_tot, YZ_Matrices.Znod, Mesh_Matrices.Gzcomplex);
    otherwise
        disp([mfilename, ': State Estimation Type ', SE_type, ' not available']);
end
[V, I] = Available_VI(V, I, Network_param.topology);

end

function [V, I] = Available_VI(V, I, topology)
    bNodes = topology.node_3ph(:) > 0;
    bBranches = topology.branch_3ph(:) > 0;

    V.Real(~bNodes) = 0;
    V.Imag(~bNodes) = 0;

    I.Real(~bBranches) = 0;
    I.Imag(~bBranches) = 0;
end

    