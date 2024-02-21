function[parameters, pseudomeas] = Network_parameters_to_SE(Network, max_nodes)
% The routine returns the network parameters and the pseudo-measurements
% Input
% -Network it is the struct defined in a file in ./Networks
%  max_nodes, if the Netwrok can be reduced, max_nodes is the maximum
%  number of nodes that the reduced network has
% Output
% -parameters
%  parameters.topology 
%  parameters.Zbr_set
%  parameters.nominal_linedata
%  parameters.mesh
% -pseudomeas
%  parameters.P
%  parameters.Q

if ~exist('max_nodes') || isempty(max_nodes) || max_nodes == -1
    max_nodes = Network.nodes.num_nodes;
end

%Base values, like base apparent power (VA) and base voltage (V).
BaseVA = Network.base.BaseVA; % kVA
BaseV = Network.base.BaseV;  % kV

avaible_nodes_reduction = Network.nodes.Available_nodes_reduction;
num_nodes = max(avaible_nodes_reduction(avaible_nodes_reduction <= max_nodes));
if isempty(num_nodes)
    num_nodes = min(avaible_nodes_reduction);
end

num_branches = min (sum(Network.branches.From <= num_nodes), sum(Network.branches.To <= num_nodes));

branch_first = Network.branches.From(1 : num_branches);
branch_end   = Network.branches.To(1 : num_branches);
if isfield(Network.branches, 'branch_cod')
    branch_cod = Network.branches.branch_cod(1: num_branches);
else
    branch_cod = [1: num_branches]';
end

parameters.topology.num_nodes = num_nodes;
parameters.topology.num_branches = num_branches;

bMesh = num_branches ~= num_nodes - 1;
if bMesh
    if length(unique([Network.branches.From(1 : num_nodes - 1), Network.branches.To(1 : num_nodes - 1)])) ~= num_nodes
        disp([mfilename, ': it is assumed that meshes branches have been added upon a radial network, please reorganize the network description']);
        assert(0);
    end  
end


branch_first_A = zeros(size(branch_first));
branch_first_B = zeros(size(branch_first));
branch_first_C = zeros(size(branch_first));
branch_end_A = zeros(size(branch_end));
branch_end_B = zeros(size(branch_end));
branch_end_C = zeros(size(branch_end));
branch_cod_A = zeros(num_branches,1);
branch_cod_B = branch_cod_A;
branch_cod_C = branch_cod_A;

branch_length = Network.branches.branch_length(1 : num_branches);       

branch_config = Network.branches.branch_config(1 : num_branches);


% il primo ramo è sempre trifase  = 1;
branch_idx_A = 1;
branch_idx_B = 1;
branch_idx_C = 1;
branch_cod_A(1) = 1;
branch_cod_B(1) = 1;
branch_cod_C(1) = 1;
isA = ones(size(branch_cod_A));
isB = ones(size(branch_cod_B));
isC = ones(size(branch_cod_C));
parameters.Zbr_set = cell(num_branches, 1);
Zbranches = zeros(num_branches, 3);
[Z3, ~] = branch_impedance(branch_impedance_per_length(branch_config(1)), branch_length(1));
% parameters.Zbr_set{1} = Z3 * (BaseVA / 1000) / (BaseV^2);
% Zbranches(1, :) = diag(Z3);
Z3 = Z3 * (BaseVA / 1000) / (BaseV^2);
parameters.Zbr_set{1} = Z3;
Zbranches(1, :) = diag(Z3);
for br_idx = 2 : num_branches
    
    [Z3, isPhs] = branch_impedance(branch_impedance_per_length(branch_config(br_idx)), branch_length(br_idx));
    isA(br_idx) = isPhs(1);
    isB(br_idx) = isPhs(2);
    isC(br_idx) = isPhs(3);
    
    Z3 = Z3 * (BaseVA / 1000) / (BaseV ^ 2);
    if (isA(br_idx))
        branch_cod_A(br_idx) = branch_idx_A + 1;
        branch_idx_A = branch_idx_A + 1;
    end
    if (isB(br_idx))
        branch_cod_B(br_idx) = branch_idx_B + 1;
        branch_idx_B = branch_idx_B + 1;
    end
    if (isC(br_idx))
        branch_cod_C(br_idx) = branch_idx_C + 1;
        branch_idx_C = branch_idx_C + 1;
    end
    Zbranches(br_idx, :) = diag(Z3);
    parameters.Zbr_set{br_idx} = Z3;
end

% kW e kVAR

P1 = Network.load.P(:, 1);
P2 = Network.load.P(:, 2);
P3 = Network.load.P(:, 3);
Q1 = Network.load.Q(:, 1);
Q2 = Network.load.Q(:, 2);
Q3 = Network.load.Q(:, 3);

P1_red = P1(1 : num_nodes - 1);
P1_red(num_nodes - 1) = sum(P1(num_nodes - 1: end));
P2_red = P2(1 : num_nodes - 1);
P2_red(num_nodes - 1) = sum(P2(num_nodes - 1: end));
P3_red = P3(1 : num_nodes - 1);
P3_red(num_nodes - 1) = sum(P3(num_nodes - 1: end));
Q1_red = Q1(1 : num_nodes - 1);
Q1_red(num_nodes - 1) = sum(Q1(num_nodes - 1: end));
Q2_red = Q2(1 : num_nodes - 1);
Q2_red(num_nodes - 1) = sum(Q2(num_nodes - 1: end));
Q3_red = Q3(1 : num_nodes - 1);
Q3_red(num_nodes - 1) = sum(Q3(num_nodes - 1: end));

P = [P1_red, P2_red, P3_red];
Q = [Q1_red, Q2_red, Q3_red];

% p.u. 5 MVA
P = P / BaseVA;
Q = Q / BaseVA;
P = -P; % iniettate
Q = -Q; % iniettate

A_node_indices = zeros(num_nodes, 1);
B_node_indices = zeros(num_nodes, 1);
C_node_indices = zeros(num_nodes, 1);
% First node is always 3Ph
A_node_indices(1) = 1;
B_node_indices(1) = 1;
C_node_indices(1) = 1;
A_idx = 1;
B_idx = 1;
C_idx = 1;
A_node_indices_reverse = [1];
B_node_indices_reverse = [1];
C_node_indices_reverse = [1];
for node_idx = 2:num_nodes
    Ato = isA(branch_end == node_idx);        % phase A presence vector of branches that arrive to node node_idx (just 1 is needed)
    Bto = isB(branch_end == node_idx);        % phase A presence vector of branches that arrive to node node_idx (just 1 is needed)
    Cto = isC(branch_end == node_idx);        % phase A presence vector of branches that arrive to node node_idx (just 1 is needed)
    if (Ato(1))                               % a A-Ph branch arrives 
       A_idx = A_idx + 1;
       A_node_indices(node_idx) = A_idx;
       A_node_indices_reverse = [A_node_indices_reverse; node_idx];
    end
    if (Bto(1)) % a A-Ph branch arrives 
       B_idx = B_idx + 1;
       B_node_indices(node_idx) = B_idx;
       B_node_indices_reverse = [B_node_indices_reverse; node_idx];
    end
    if (Cto(1)) % a A-Ph branch arrives 
       C_idx = C_idx + 1;
       C_node_indices(node_idx) = C_idx;
       C_node_indices_reverse = [C_node_indices_reverse; node_idx];
    end
end

parameters.topology.node_3ph = [A_node_indices, B_node_indices, C_node_indices];
num_A_nodes = length(A_node_indices_reverse);
num_B_nodes = length(B_node_indices_reverse);
num_C_nodes = length(C_node_indices_reverse);
parameters.topology.node_3ph_reverse = zeros(max([num_A_nodes,num_B_nodes,num_C_nodes]),3);
parameters.topology.node_3ph_reverse(1:num_A_nodes,1) = A_node_indices_reverse;
parameters.topology.node_3ph_reverse(1:num_B_nodes,2) = B_node_indices_reverse;
parameters.topology.node_3ph_reverse(1:num_C_nodes,3) = C_node_indices_reverse;

branch_first_A(1) = 1;
branch_first_B(1) = 1;
branch_first_C(1) = 1;
for br_idx = 1:num_branches
    from = branch_first(br_idx);
    to = branch_end(br_idx);
    if (isA(br_idx))
        branch_first_A(br_idx) = A_node_indices(from);
        branch_end_A(br_idx) = A_node_indices(to);
    end
    if (isB(br_idx))
        branch_first_B(br_idx) = B_node_indices(from);
        branch_end_B(br_idx) = B_node_indices(to);
    end
    if (isC(br_idx))
        branch_first_C(br_idx) = C_node_indices(from);
        branch_end_C(br_idx) = C_node_indices(to);
    end
end


parameters.topology.num_state_nodes = num_A_nodes + num_B_nodes + num_C_nodes;
parameters.topology.branch_3ph = [isA, isB, isC];

parameters.nominal_linedata = table(branch_first, branch_end, ... % 1, 2
    branch_first_A, branch_first_B, branch_first_C, ... % 3, 4, 5
    branch_end_A, branch_end_B, branch_end_C, ... % 6, 7, 8
    branch_cod, branch_cod_A, branch_cod_B, branch_cod_C, ... % 9, 10, 11, 12
    isA, isB, isC, branch_config); % 13, 14, 15, 16, 17, 18, 19
pseudomeas.P = P;
pseudomeas.Q = Q;

 %%% Generation eq_curr_matrix

PA = P(isA(1:num_nodes-1)==1, 1);   QA = Q(isA(1:num_nodes-1)==1, 1);
PB = P(isB(1:num_nodes-1)==1, 2);   QB = Q(isB(1:num_nodes-1)==1, 2);
PC = P(isC(1:num_nodes-1)==1, 3);   QC = Q(isC(1:num_nodes-1)==1, 3);

branch_first_A = branch_first_A(isA==1);   branch_end_A = branch_end_A(isA==1);  
branch_first_B = branch_first_B(isB==1);   branch_end_B = branch_end_B(isB==1);
branch_first_C = branch_first_C(isC==1);   branch_end_C = branch_end_C(isC==1);

num_branchesA = sum(branch_cod_A > 0);
num_branchesB = sum(branch_cod_B > 0);
num_branchesC = sum(branch_cod_C > 0);

zero_injA = find(PA==0 & QA==0);
zero_injB = find(PB==0 & QB==0);
zero_injC = find(PC==0 & QC==0);

zero_inj_revA = flipud(zero_injA);
zero_inj_revB = flipud(zero_injB);
zero_inj_revC = flipud(zero_injC);

matrixA = eye(num_branchesA,num_branchesA);
matrixB = eye(num_branchesB,num_branchesB);
matrixC = eye(num_branchesC,num_branchesC);

for i=1:length(zero_inj_revA)
    cod = find(branch_first_A == zero_inj_revA(i)+1);
    for j=1:length(cod)
        matrixA(zero_inj_revA(i),:) = matrixA(zero_inj_revA(i),:) + matrixA(cod(j,1),:);
    end
    cod = find(branch_end_A == zero_inj_revA(i)+1);
    for j=1:length(cod)
        if cod(j) ~= zero_inj_revA(i)
           matrixA(zero_inj_revA(i),:) = matrixA(zero_inj_revA(i),:) - matrixA(cod(j,1),:);
        end
    end
end
for i=1:length(zero_inj_revB)
    cod = find(branch_first_B == zero_inj_revB(i)+1);
    for j=1:length(cod)
        matrixB(zero_inj_revB(i),:) = matrixB(zero_inj_revB(i),:) + matrixB(cod(j,1),:);
    end
    cod = find(branch_end_B == zero_inj_revB(i)+1);
    for j=1:length(cod)
        if cod(j) ~= zero_inj_revB(i)
           matrixB(zero_inj_revB(i),:) = matrixB(zero_inj_revB(i),:) - matrixB(cod(j,1),:);
        end
    end
end
for i=1:length(zero_inj_revC)
    cod = find(branch_first_C == zero_inj_revC(i)+1);
    for j=1:length(cod)
        matrixC(zero_inj_revC(i),:) = matrixC(zero_inj_revC(i),:) + matrixC(cod(j,1),:);
    end
    cod = find(branch_end_C == zero_inj_revC(i)+1);
    for j=1:length(cod)
        if cod(j) ~= zero_inj_revC(i)
           matrixC(zero_inj_revC(i),:) = matrixC(zero_inj_revC(i),:) - matrixC(cod(j,1),:);
        end
    end
end

matrixA(:,zero_inj_revA) = [];
matrixB(:,zero_inj_revB) = [];
matrixC(:,zero_inj_revC) = [];
dimA = size(matrixA);
dimB = size(matrixB);
dimC = size(matrixC);
matrix2A = zeros(dimA);
matrix2B = zeros(dimB);
matrix2C = zeros(dimC);

matrixAB = zeros(size(matrixA,1),size(matrixB,2));
matrixAC = zeros(size(matrixA,1),size(matrixC,2));
matrixBA = zeros(size(matrixB,1),size(matrixA,2));
matrixBC = zeros(size(matrixB,1),size(matrixC,2));
matrixCA = zeros(size(matrixC,1),size(matrixA,2));
matrixCB = zeros(size(matrixC,1),size(matrixB,2));

parameters.mesh.eq_curr_matrix = eye(2*(dimA+dimB+dimC+3));
parameters.mesh.eq_curr_matrix(7:end,7:end) = [matrixA, matrixAB, matrixAC, matrix2A, matrixAB, matrixAC;
                               matrixBA, matrixB, matrixBC, matrixBA, matrix2B, matrixBC;
                               matrixCA, matrixCB, matrixC, matrixCA, matrixCB, matrix2C;
                               matrix2A, matrixAB, matrixAC, matrixA, matrixAB, matrixAC;
                               matrixBA, matrix2B, matrixBC, matrixBA, matrixB, matrixBC;
                               matrixCA, matrixCB, matrix2C, matrixCA, matrixCB, matrixC];
                           
parameters.mesh.zero_inj_cell = cell(1,3);
parameters.mesh.zero_inj_cell{1,1} = zero_injA;
parameters.mesh.zero_inj_cell{1,2} = zero_injB;
parameters.mesh.zero_inj_cell{1,3} = zero_injC;

end

function [Z3, bPh] = branch_impedance(Z3_per_unit, branch_length)
    bPh = sum(abs(Z3_per_unit), 2) > 0;
    Z3 = branch_length * (Z3_per_unit + triu(Z3_per_unit, 1).');
end
