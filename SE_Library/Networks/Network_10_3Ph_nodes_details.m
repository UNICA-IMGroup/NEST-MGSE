function Network = Network_10_3Ph_nodes_details()
%%
%Network Topology
Network.nodes.num_nodes = 10;
Network.nodes.Available_nodes_reduction = [10];  

Network.branches.From = [1, 2, 3, 2, 2, 6, 7, 6, 7]';
Network.branches.To   = [2, 3, 4, 5, 6, 7, 8, 9, 10]';

%%
%Network Parameters
%branch length in miles
Network.branches.branch_length = [2000, 500, 300, 500, 2000, 300, 300, 500, 800]/5280; % miles (5280 feet per mile)
             
%branch cofig, the network is composed by different physical connectors
%whose physical parameters are detailed in the file
%branch_impedence_from_config.com
%config: code of the type of branch
Network.branches.branch_config = [601, 603, 603, 602, 601, 604, 605, 606, 607]';
     
%%
%Network base values
Network.base.BaseVA = 5000; % kVA
Network.base.BaseV = 4.16;  % kV
% BaseV = 4.16/sqrt(3);  % kV

%Active Power
Network.load.P = [8.5, 33, 58.5; 0, 170, 170; 0, 230, 230; 160, 120, 120; 385 + 8.5 + 0, 385 + 33 + 0, 385 + 58.5 + 170; 0, 0, 0; 0, 0, 170;       485,       68,       290;       128, 0, 0];

%Reactive Power
Network.load.Q = [5,   19, 34;   0, 125, 125; 0, 132, 132; 110, 90,  90;  220 + 5   + 0, 220 + 19 + 0, 220 + 34   + 151; 0, 0, 0; 0, 0, 80 - 200;  190 - 200, 60 - 200, 212 - 200; 86,  0, 0];
