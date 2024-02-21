function [YZ_Matrices, Mesh_Matrices] = Impedence_Matrices_to_SE(Network_param)

nominal_linedata = table2array(Network_param.nominal_linedata);

[YZ_Matrices.Y_state, YZ_Matrices.Y_tot, YZ_Matrices.Znod, YZ_Matrices.Znod_state] ...
    = Ymatrix_3Ph_to_SE(nominal_linedata, Network_param.Zbr_set, Network_param.topology.num_nodes);    %%% calcolo matrice ammettenza della rete

[Mesh_Matrices.Gz, Mesh_Matrices.Gz_state, Mesh_Matrices.Gzmod, Mesh_Matrices.Gzangle, Mesh_Matrices.Gzmod_state, Mesh_Matrices.Gzangle_state, Mesh_Matrices.Gzcomplex] ...
    = Mesh_matrix_3Ph(Network_param.topology.num_nodes, nominal_linedata, Network_param.topology.num_branches, Network_param.Zbr_set, ...
    Network_param.mesh.zero_inj_cell, Network_param.mesh.eq_curr_matrix);      %%% calcolo delle matrici contenenti le impedenze delle maglie



