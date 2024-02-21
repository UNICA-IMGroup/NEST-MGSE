%Availables Network description files
% Network_118_3Ph_nodes_details;
% Network_118_3Ph_nodes_meshed_details;
% Network_10_3Ph_nodes_details;
% Network_10_3Ph_nodes_meshed_details
Load_NetworkDataFunction = @Network_10_3Ph_nodes_details;

% the network Network_118_3Ph_nodes_details has 118 nodes but it is
% possibile to use a reduced version of the network
%Network.nodes.Available_nodes_reduction = [8, 15, 54 : 56, 59, 62, 69, 99, 104, 108, 111, 113, 114, 116, 118];  
max_nodes = [];


Load_MeasurementConfiguration = @MeasurementConfiguration_to_SE;

% Measurement placement configuration - file MeasurementConfiguration_to_SE.m 
% 'basic placement' -> only the first node and the first branch are
% monitored
% '3 nodes placement', in addiction also nodes 7 and 9 are monitored.
% it is possible to add placement_config specifications in the file MeasurementConfiguration_to_SE.m 
%placement_config_type = {'default', '3 nodes placement' };
placement_config = 'default';

%Measurement Uncertainty config - file Measurement_uncertainty.m 
% unc_config_type = {'nominal', 'default'}
unc_config = 'nominal';

%Distributed Generation configuration. Each node connected to DG is
%monitored with a Voltage and Power injection measurements
% It is possible to choose among:
% 'no DG' no DG is avaible
% 'last node' the DG is placed in the last node
% 'specific' it is possible to choose a list of nodes with DG generation
%DG_placement_type = {'no DG', 'last node', 'specific'};
DG_placement = 'last node';

if strcmp(DG_placement, 'specific')
    DG_position = [4, 6];
end

%number of Monte Carlo trials
MC_iterations = 1000;

% each MC trial corresponds to a different operative case of the network
% obtained by extracting PQ from a normal distrition centered in nominal PQ
% values (definded in the actual Load_NetworkDataFunction)
% having standard deviation deltaPQ / 3.
deltaPQ = 50; 

%plots

%if true, the Network graph is ploted
bGraphPlot = true;

%at the end of the execution, the Voltage profiles, relative to a specific MC iteration
%(reference and %estimated) is plotted 
mciter_to_plot = 3;