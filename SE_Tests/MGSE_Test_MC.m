function MGSE_Test_MC(fileTestConfig)

close all;
clear, clc;
bSaveData = true;

rng(0);
clear Newton_Raphson_Power_Flow_to_SE
clear State_Estimation
   
bPlotStatistics = true;

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


%% Montecarlo Simulation
% Type of State Estimator: Branch currents
SE_type = ["BC rect"];

% Results initialization
estimates_MC = estimatesInitialize(num_nodes, num_branches, 3,  MC_iterations); 
reference_values_MC = VIS_inizialize_MC(num_nodes, num_branches, 3, MC_iterations);
bInitialize = true;

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

    if bInitialize
        V0 = ones (3 * num_nodes, 1);
    end
    
        tic 
    %% State Estimation
    [V, I] = MG_State_Estimation(Measurements, pseudomeas, Measurements_param, Network_param, bInitialize, V0);
    V0 = V.Module(:);
    bInitialize = false;
        
    %% results
    estimates_MC.Time(mciter) = toc; 
    estimates_MC.V.Module(:, mciter) = V.Module;
    estimates_MC.V.Angle(:, mciter) = V.Angle;
    estimates_MC.I.Module(:, mciter) = I.Module;
    estimates_MC.I.Angle(:, mciter) = I.Angle;
    estimates_MC.V.Unc.Module(:, mciter) = V.Unc.Module;
    estimates_MC.V.Unc.Angle(:, mciter) = V.Unc.Angle;
    estimates_MC.I.Unc.Module(:, mciter) = I.Unc.Module;
    estimates_MC.I.Unc.Angle(:, mciter) = I.Unc.Angle;

    if(mciter == mciter_to_plot)
        estimates_to_Plot.V = V;
        estimates_to_Plot.I = I;
        reference_values_to_Plot = reference_values;
        general_title = strcat('#Nodes:', num2str(num_nodes), ' - #Branches:', num2str(num_branches), ...
            ' - #Meas Points: ', num2str(nMis), ' - #Pseudomeasurements Uncertainty:', num2str(deltaPQ), ' %');
    end

end

figures(iFigure) = plotElapsedTime(estimates_MC.Time, SE_type);
iFigure = iFigure + 1;

hV = plotVoltages_3Ph_Advanced(reference_values_MC.V, estimates_MC.V, Network_param.topology, MC_iterations, mciter_to_plot, bPlotStatistics);
nF = length(hV);
figures(iFigure : iFigure + nF - 1) = hV;
iFigure = iFigure + nF;

hC = plotCurrents_3Ph_Advanced(reference_values_MC.I, estimates_MC.I, Network_param.topology, MC_iterations, mciter_to_plot, bPlotStatistics);
nF = length(hC);
figures(iFigure : iFigure + nF - 1) = hC;
iFigure = iFigure + nF;

hVI = plotVI_3Ph_Advanced(estimates_to_Plot.V, estimates_to_Plot.I, ...
    reference_values_to_Plot.V(:), reference_values_to_Plot.I(:), Network_param.topology, general_title);
figures(iFigure) = hVI; 

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

end

function estimates = estimatesInitialize(num_nodes, num_branches, num_phases, MC)

estimates.V.Module = zeros(num_nodes * num_phases, MC);
estimates.V.Angle = zeros(num_nodes * num_phases, MC);
estimates.I.Module = zeros(num_branches * num_phases, MC);
estimates.I.Angle = zeros(num_branches * num_phases, MC);

estimates.V.Unc.Module = zeros(num_nodes * num_phases, MC);
estimates.V.Unc.Angle = zeros(num_nodes * num_phases, MC);
estimates.I.Unc.Module = zeros(num_branches * num_phases, MC);
estimates.I.Unc.Angle = zeros(num_branches * num_phases, MC);

estimates_MC.Time = zeros(MC, 1);
end
