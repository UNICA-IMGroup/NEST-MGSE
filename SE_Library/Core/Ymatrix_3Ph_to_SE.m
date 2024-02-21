function[Y_state, Y_tot, Znod, Znod_state] = Ymatrix_3Ph_to_SE(nominal_linedata, Zbr_set, num_nodes)
% function[Y_state] = Ymatrix_3Ph(num_nodes, nominal_linedata, Zbr_set)
% Ammettance matrix for a three phase network. The matrix is 
% [N X M] where N is the number of existing branch-current variables (they
% are not 3*(num_nodes-1) - in a radial system - because the system is not 3Ph in every branch)
% and M is the number of existing voltage variables (they
% are not 3*(num_nodes) because the system is not 3Ph in every node).
% in: 
% - num_nodes
% - nominal_linedata: all the needed topology data
% - Zbr_set: set that contains a 3x3 branch impedance for each branch
% out:
% - Y_state: ammettance matrix

num_branches = size(nominal_linedata, 1);
% nominal_linedata.branches = nominal_linedata(:, [1 2]);
% nominal_linedata
branch_first = nominal_linedata(:,1);
branch_end = nominal_linedata(:,2);
branch_first_A = nominal_linedata(:,3);
branch_first_B = nominal_linedata(:,4);
branch_first_C = nominal_linedata(:,5);
branch_end_A = nominal_linedata(:,6);
branch_end_B = nominal_linedata(:,7);
branch_end_C = nominal_linedata(:,8);

isA = nominal_linedata(:,13);
isB = nominal_linedata(:,14);
isC = nominal_linedata(:,15);
num_nodes_A = max(branch_end_A);
num_nodes_B = max(branch_end_B);
num_nodes_C = max(branch_end_C);
% num_br_A = max(branch_cod_A);
% num_br_B = max(branch_cod_B);
% num_br_C = max(branch_cod_C);

Y_state = zeros(num_nodes_A + num_nodes_B + num_nodes_C, num_nodes_A + num_nodes_B + num_nodes_C);
Y_tot = zeros(3*num_nodes,3*num_nodes);
for L=1:num_branches              %%% calcolo componenti fuori dalla diagonale della matrice di ammettenza
    Z3true = Zbr_set{L};
    remove_index = 3;
    if (~isA(L))
        Z3true(remove_index - 2, :) = [];
        Z3true(:, remove_index - 2) = [];
        remove_index = remove_index - 1;
    end
    if (~isB(L))
        Z3true(remove_index - 1, :) = [];
        Z3true(:, remove_index - 1) = [];
        remove_index = remove_index - 1;
    end
    if (~isC(L))
        Z3true(remove_index, :) = [];
        Z3true(:, remove_index) = [];
    end
    m = branch_first(L);
    n = branch_end(L);
    mm = [m; m+num_nodes; m+2*num_nodes];
    nn = [n; n+num_nodes; n+2*num_nodes];
    iA = branch_first_A(L);
    jA = branch_end_A(L);
    iB = branch_first_B(L);
    jB = branch_end_B(L);
    iC = branch_first_C(L);
    jC = branch_end_C(L);
    
    Y3 = inv(Z3true);
    Y3tot = Y3;
    if (~isA(L))
        dim = size(Y3tot,1);
        Y3tot = [zeros(dim,1),Y3tot];
        Y3tot = [zeros(1,dim+1);Y3tot];
    end
    if (~isB(L))
        dim = size(Y3tot,1);
        if dim == 1
            Y3tot = [Y3tot,zeros(dim,1)];
            Y3tot = [Y3tot;zeros(1,dim+1)];
        else
            Y3tot = [Y3tot(:,1),zeros(dim,1),Y3tot(:,2)];
            Y3tot = [Y3tot(1,:);zeros(1,dim+1);Y3tot(2,:)];
        end
    end
    if (~isC(L))
        Y3tot = [Y3tot,zeros(2,1)];
        Y3tot = [Y3tot;zeros(1,3)];
    end
       
    node_phase_i = [];
    node_phase_j = [];
    if (iA ~= 0 && jA ~= 0)
        node_phase_i = [node_phase_i, iA];
        node_phase_j = [node_phase_j, jA];
    end
    if (iB ~= 0 && jB ~= 0)
        node_phase_i = [node_phase_i, iB + num_nodes_A];
        node_phase_j = [node_phase_j, jB + num_nodes_A];
    end
    if (iC ~= 0 && jC ~= 0)
        node_phase_i = [node_phase_i, iC + num_nodes_A + num_nodes_B];
        node_phase_j = [node_phase_j, jC + num_nodes_A + num_nodes_B];
    end
    
    if L == 68
        a = L;
    end
    
    Y_state(node_phase_i, node_phase_j) = -Y3;
    Y_state(node_phase_j, node_phase_i) = -Y3; % se traspongo è lo stesso
%     Y_state(j,  i)    = Y_state(i,j);
    Y_state(node_phase_i, node_phase_i) = Y_state(node_phase_i, node_phase_i) + Y3;
    Y_state(node_phase_j, node_phase_j) = Y_state(node_phase_j, node_phase_j) + Y3;
    Y_tot(mm,nn) = - Y3tot;
    Y_tot(nn,mm) = - Y3tot;
    Y_tot(mm,mm) = Y_tot(mm,mm) + Y3tot;
    Y_tot(nn,nn) = Y_tot(nn,nn) + Y3tot;
end
Znod = zeros(3*(num_nodes-1), 3*num_branches);       %%%    matrice delle impedenze da usare con le misure di tensione Vmag_nod  

% solo i primi num_nodes - 1 branches cioè quelli dell'albero
for i=1:num_nodes-1
    a=branch_first(i)-1; % a rappresenta il precedente nodo nella lista
    if a==0 % era il primo ramo che parte da 1 --> 
        Znod([i, i + num_nodes - 1, i + 2*(num_nodes - 1)],:)=Znod([i, i + num_nodes - 1, i + 2*(num_nodes - 1)], :);
    else
        Znod([i, i + num_nodes - 1, i + 2*(num_nodes - 1)],:)=Znod([a, a + num_nodes - 1, a + 2*(num_nodes - 1)],:);
    end
    Z3 = Zbr_set{i};
    Znod([i, i + num_nodes - 1, i + 2*(num_nodes - 1)], [i, i + num_branches, i + 2*num_branches]) = Z3;
end
isAA = isA>0;
isAA1 = isAA(1:num_nodes-1);
isBB = isB>0;
isBB1 = isBB(1:num_nodes-1);
isCC = isC>0;
isCC1 = isCC(1:num_nodes-1);
isABC = [isAA;isBB;isCC];
isABC1 = [isAA1;isBB1;isCC1];
Znod_state = Znod(isABC1,isABC);
