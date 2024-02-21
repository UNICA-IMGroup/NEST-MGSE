function [reference_values] = Newton_Raphson_Power_Flow_to_SE(Network_param, pseudomeas, V_init)

persistent Y_state;

V_init.V_abs = ones (3 * Network_param.topology.num_nodes, 1);                                                           %%% inizializzazione dei vettori delle tensioni ai nodi
V_init.V_theta = [zeros(Network_param.topology.num_nodes, 1); 
    ones(Network_param.topology.num_nodes, 1) * ( 2 * pi / 3); 
    ones(Network_param.topology.num_nodes, 1) * (-2 * pi / 3)];

if isempty(Y_state)
    [Y_state] = Ymatrix_3Ph_to_SE(table2array(Network_param.nominal_linedata), Network_param.Zbr_set, Network_param.topology.num_nodes);    %%% calcolo matrice ammettenza della rete
end

[reference_values.V, reference_values.I, reference_values.S_br, reference_values.S_inj, ~] ...
    = Power_flow_New_Raph_3Ph(Network_param.topology.num_nodes, Network_param.topology.node_3ph, Network_param.topology.node_3ph_reverse, ...
                Y_state, pseudomeas, table2array(Network_param.nominal_linedata), Network_param.Zbr_set, 1.01 * V_init.V_abs, V_init.V_theta);     
