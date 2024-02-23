
function SE_Test_MC(fileTestConfig)

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

%% Graphic rappresentation of network and measurement system
if bGraphPlot    
    figures(1) = PlotNetwork_to_SE(Network_param.nominal_linedata, Measurements_param.meas_indices);
end

%% Montecarlo Simulation
% Types of State Estimator, Branch currents and Node Voltages
SE_type_Group = ["BC rect", "NV rect"];

% Results initialization
estimates_MC = estimatesInitialize(num_nodes, num_branches, 3, length(SE_type_Group), MC_iterations); 
reference_values_MC = VIS_inizialize_MC(num_nodes, num_branches, 3, MC_iterations);

%each MC_iterations trial corresponds to a different operative case of the network
for mciter = 1:MC_iterations                                                             %%% inizio delle iterazioni MonteCarlo                   
    if mod(mciter, 100)==0                                                     
        mciter
    end
        
    % PQ extraction from a normal distrition centered in nominal P values
    % having standard deviation deltaPQ / 3.
    pseudomeas_trial = [pseudomeas.P, pseudomeas.Q] .* (1 + deltaPQ / (3 * 100) .*(1 + randn(size([pseudomeas.P, pseudomeas.Q]))));

    %The power flow compute reference values for V, I, S_br, S_inj
    reference_values = Newton_Raphson_Power_Flow_to_SE(Network_param, pseudomeas_trial); 

    reference_values_MC.V(:,  mciter) = reference_values.V(:); reference_values_MC.I(:,  mciter) = reference_values.I(:); 
    reference_values_MC.S_br(:,  mciter) = reference_values.S_br(:); reference_values_MC.S_inj(:,  mciter) = reference_values.S_inj(:);
    
    %% Simulation of Measurements: reference values are transformed in
    % measurement using Measurement parameters 
    Measurements = AddMeasurementErrors_SE(reference_values, Measurements_param);
   
    iter_method = 1;
    bInitialize = mciter == 1;
    for SE_type = SE_type_Group
        tic 
        %% State Estimation
        [V, I] = State_Estimation(Measurements, pseudomeas, Measurements_param, Network_param, SE_type, bInitialize);
        
        %% results
        estimates_MC.Time(iter_method, mciter) = toc; 
        estimates_MC.V(:, iter_method, mciter) = complex(V.Real, V.Imag);
        estimates_MC.I(:, iter_method, mciter) = complex(I.Real, I.Imag);
        iter_method = iter_method + 1;
        %%
    end
end

figures(2) = plotElapsedTime(estimates_MC.Time, SE_type_Group);
figures([3 : 6]) = plotVoltages_3Ph(reference_values_MC.V, estimates_MC.V, Network_param.topology, SE_type_Group, MC_iterations, mciter_to_plot);
figures([7 : 10]) = plotCurrents_3Ph(reference_values_MC.I, estimates_MC.I, Network_param.topology, SE_type_Group, MC_iterations, mciter_to_plot);

Rem_Paths(additional_Paths);

if bSaveData
    final_results_folder = strcat(results_folder, '/', erase(fileTestConfig, '.m'), '_', datestr(now,'mm-dd_HH_MM_SS'));
    mkdir(final_results_folder);
    filename = strcat(erase(fileTestConfig, '.m'), '_results');
    savefig(figures, strcat(final_results_folder, '/', filename));
    clear figures;
    save(strcat(final_results_folder, '/', filename, '.mat'));
    
    copyfile(strcat(mfilename('fullpath'),'.m'), final_results_folder);
end
