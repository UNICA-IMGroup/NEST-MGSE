function MGSE_Test(fileTestConfig)

close all;
clear, clc;
bSaveData = true;

rng(0);
clear Newton_Raphson_Power_Flow_to_SE
clear State_Estimation

if ~exist('fileTestConfig', 'var')
     fileTestConfig = 'Test_MC_config_3';
else
    fileTestConfig = erase(fileTestConfig, '.m');
end

results_folder = './results';
additional_Paths = Add_Paths;
eval(fileTestConfig)


%% Load network data
Network = Load_NetworkDataFunction();
[Network_param, pseudomeas] ...
    = Network_parameters_to_SE(Network, max_nodes);
num_nodes = Network_param.topology.num_nodes;
num_branches = Network_param.topology.num_branches;

%% Load measurement system
switch DG_placement
    case 'last node'        
        DG_position = num_nodes;
    case 'specific'
        DG_position = DG_position(DG_position <= num_nodes);
    otherwise
        DG_position = [];
end
Measurements_param = Load_MeasurementConfiguration(Network_param.topology.num_nodes, Network_param.topology.num_branches, DG_position, placement_config, unc_config);

iFigure = 1;
%% Graphic rappresentation of network and measurement system
if bGraphPlot    
    [figures(iFigure), nMis] = PlotNetwork_to_SE(Network_param.nominal_linedata, Measurements_param.meas_indices);
    iFigure = iFigure + 1;
else
     monitoredNodes = unique([Measurements_param.meas_indices.PQ_inj_idx; Measurements_param.meas_indices.Vmag_nod_idx]);
     monitoredBranches_index = unique([Measurements_param.meas_indices.Imag_br_idx; Measurements_param.meas_indices.PQ_br_idx]);
     nMis = length(monitoredNodes) + length(monitoredBranches_index);
end

general_title = strcat('#Nodes:', num2str(num_nodes), ' - #Branches:', num2str(num_branches), ...
            ' - #Meas Points: ', num2str(nMis), ' - #Pseudomeasurements Uncertainty:', num2str(deltaPQ), ' %');

%% Montecarlo Simulation
MC_iterations = 10;
% Results initialization

bInitialize = true;

%each MC_iterations trial corresponds to a different operative case of the network
for mciter = 1:MC_iterations                                                             %%% inizio delle iterazioni MonteCarlo                   
        
    % PQ extraction from a normal distrition centered in nominal P values
    % having standard deviation deltaPQ / 3.
    pseudomeas_trial = [pseudomeas.P, pseudomeas.Q] .* (1 + deltaPQ / (3 * 100) .*(1 + randn(size([pseudomeas.P, pseudomeas.Q]))));

    %The power flow compute reference values for V, I, S_br, S_inj
    reference_values = Newton_Raphson_Power_Flow_to_SE(Network_param, pseudomeas_trial); 
    
    %% Simulation of Measurements: reference values are transformed in
    % measurement using Measurement parameters 
    Measurements = AddMeasurementErrors_SE(reference_values, Measurements_param);

    if bInitialize
        V0 = ones (3 * num_nodes, 1);
    end
    
        tic 
    %% State Estimation
    [V, I] = MG_State_Estimation(Measurements, pseudomeas, Measurements_param, Network_param, bInitialize, V0);

    V0 = V.Module(:);

    bInitialize = false;

    hVI = plot_MGSE(V, I, Network_param.topology, general_title);   

    figures(iFigure) = hVI;
    iFigure = iFigure + 1;
end
end

