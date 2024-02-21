function[zdata] = Zdata_3Ph_to_SE(nominal_linedata_, pseudomeas, values, meas, bUniform)

if isstruct(nominal_linedata_)
    nominal_linedata = struct2array(nominal_linedata_);
else
    nominal_linedata = nominal_linedata_;
end

if isstruct(pseudomeas)
    P = pseudomeas.P;
    Q = pseudomeas.Q;
else
    P = pseudomeas(:,1:3);
    Q = pseudomeas(:,4:6);
end
    
if ~exist('bUniform', 'var')
    bUniform = false;
end
[V, I_br, S_br, S_inj] = extract_reference_values(values);    
[Vmag_nod_idx, Imag_br_idx, Isync_magphase_br_idx, Vsync_magphase_nod_idx, PQ_br_idx, PQ_inj_idx] = extract_meas_indices(meas.meas_indices);

[std_rel_Vmag_nod, std_rel_Imag_br, std_rel_Isync_mag, std_rel_Isync_phase, std_rel_Vsync_mag, std_rel_Vsync_phase, ...
    std_rel_PQ_br, std_rel_PQ_inj, std_rel_PQ_pseudo] = extract_meas_uncertainty(meas.meas_unc, bUniform);


branch_first = nominal_linedata(:,1);
branch_end = nominal_linedata(:,2);
branch_cod = nominal_linedata(:,9);

numPi = size(P, 1);
numQi = size(Q, 1);
numPf = length(PQ_br_idx);
numQf = length(PQ_br_idx);
numV = length(Vmag_nod_idx);
numI = length(Imag_br_idx);
numVsync = length(Vsync_magphase_nod_idx);
numIsync = length(Isync_magphase_br_idx);
num_meas = numPi + numQi + numPf + numQf + numV + numI + 2 * numVsync + 2 * numIsync;

P_inj = real(S_inj(PQ_inj_idx, :));                                  %%% definizione delle misure esatte
Q_inj = imag(S_inj(PQ_inj_idx, :));
P(PQ_inj_idx - 1, :) = P_inj;
Q(PQ_inj_idx - 1, :) = Q_inj;

P_br = real(S_br(PQ_br_idx, :));                                  %%% definizione delle misure esatte
Q_br = imag(S_br(PQ_br_idx, :));
Vmag_nod = abs(V(Vmag_nod_idx, :));
Imag_br = abs(I_br(Imag_br_idx, :));
Isync_mag_br = abs(I_br(Isync_magphase_br_idx, :));
Isync_phase_br = angle(I_br(Isync_magphase_br_idx, :));
Vsync_mag_nod = abs(V(Vsync_magphase_nod_idx,:));
Vsync_phase_nod = angle(V(Vsync_magphase_nod_idx,:));

t1 = ones(numV, 1);
t2 = 2*ones(numPi, 1);
t3 = 3*ones(numQi, 1);
t4 = 4*ones(numPf, 1);
t5 = 5*ones(numQf, 1);
t6 = 6*ones(numI, 1);
t7 = 7*ones(numVsync, 1);
t8 = 8*ones(numVsync, 1);
t9 = 9*ones(numIsync, 1);
t10 = 10*ones(numIsync, 1);
type = [t1; t2; t3; t4; t5; t6; t7; t8; t9; t10];                          %%% colonna1 della matrice zdata

z_true = [Vmag_nod; P; Q; P_br; Q_br; Imag_br; Vsync_mag_nod; Vsync_phase_nod; Isync_mag_br; Isync_phase_br];   %%% colonne 2,3,4 della matrice zdata

fromPQ_inj = [2 : numPi+1]'; 
toPQ_inj = zeros(numPi,1);
brPQ_inj = zeros(numPi,1);
fromP_br = branch_first(PQ_br_idx);
toP_br = branch_end(PQ_br_idx);
brP_br = branch_cod(PQ_br_idx);
fromQ_br = branch_first(PQ_br_idx);
toQ_br = branch_end(PQ_br_idx);
brQ_br = branch_cod(PQ_br_idx);
fromV = Vmag_nod_idx;
toV = zeros(numV,1);
brV = zeros(numV,1);
fromI = branch_first(Imag_br_idx);
toI = branch_end(Imag_br_idx);
brI = branch_cod(Imag_br_idx);
fromVsync = Vsync_magphase_nod_idx;
toVsync = zeros(numVsync,1);
brVsync = zeros(numVsync,1);
fromIsync = branch_first(Isync_magphase_br_idx);
toIsync = branch_end(Isync_magphase_br_idx);
brIsync = branch_cod(Isync_magphase_br_idx);
br = [brV; brPQ_inj; brPQ_inj; brP_br; brQ_br; brI; brVsync; brVsync; brIsync; brIsync];   %%% colonna5 della matrice zdata
from = [fromV; fromPQ_inj; fromPQ_inj; fromP_br; fromQ_br; fromI; fromVsync; fromVsync; fromIsync; fromIsync];   %%% colonna6 della matrice zdata
to = [toV; toPQ_inj; toPQ_inj; toP_br; toQ_br; toI; toVsync; toVsync; toIsync; toIsync];   %%% colonna7 della matrice zdata

% PQ inj
dev_std_Pinj = std_rel_PQ_pseudo * ones(numPi,3);
dev_std_Qinj = std_rel_PQ_pseudo * ones(numQi,3);

dev_std_Pinj(PQ_inj_idx - 1, :) = std_rel_PQ_inj;
dev_std_Qinj(PQ_inj_idx - 1, :) = std_rel_PQ_inj;

dev_std_Pf = std_rel_PQ_br * ones(numPf,3);
dev_std_Qf = std_rel_PQ_br * ones(numQf,3);
dev_std_V = std_rel_Vmag_nod * ones(numV,3);
dev_std_I = std_rel_Imag_br * ones(numI,3);
dev_std_Vsync_mag = std_rel_Vsync_mag * ones(numVsync,3);
dev_std_Vsync_phase = std_rel_Vsync_phase * ones(numVsync,3);
dev_std_Isync_mag = std_rel_Isync_mag * ones(numIsync,3);
dev_std_Isync_phase = std_rel_Isync_phase *ones(numIsync,3);

%%%% fine Vincenzo 22/01/2023
dev_std_rel = [dev_std_V; dev_std_Pinj; dev_std_Qinj; dev_std_Pf; dev_std_Qf; dev_std_I;
               dev_std_Vsync_mag; dev_std_Vsync_phase; dev_std_Isync_mag; dev_std_Isync_phase];
dev_std = dev_std_rel .* abs(z_true);
abs_acc_i = (type == 8 | type == 10);
abs_acc = find(abs_acc_i);
dev_std(abs_acc,:) = dev_std_rel(abs_acc, :);                               %%% colonne 8,9,10 della matrice zdata

zdata = zeros(num_meas, 10);
zdatatrue2 = zeros(num_meas, 10);
zdatatrue2(:, 1) = type;
zdatatrue2(:, 2:4) = z_true;
zdatatrue2(:, 5) = br;
zdatatrue2(:, 6) = from;
zdatatrue2(:, 7) = to;
zdatatrue2(:, 8:10) = dev_std;

Vsyncstart = 1 + numPi + numQi + numPf + numQf + numV + numI;
Isyncstart = Vsyncstart + 2*numVsync;
zdata(1:Vsyncstart-1 , :) = zdatatrue2(1:Vsyncstart-1 , :);
for i=1:numVsync
    zdata(Vsyncstart + 2*(i-1) , :) = zdatatrue2(Vsyncstart + i-1 , :);
    zdata(Vsyncstart + 2*i - 1 , :) = zdatatrue2(Vsyncstart + i-1 + numVsync , :);
end
for i=1:numIsync
    zdata(Isyncstart + 2*(i-1) , :) = zdatatrue2(Isyncstart + i-1 , :);
    zdata(Isyncstart + 2*i - 1 , :) = zdatatrue2(Isyncstart + i-1 + numIsync , :);
end
type = zdata(:, 1);
z = zdata(:, 2:4);
br = zdata(:, 5);
from = zdata(:, 6);
to = zdata(:, 7);
dev_std = zdata(:, 8:10);

zdata = table(type, z, br, from, to, dev_std);
end

function [V, I_br, S_br, S_inj] = extract_reference_values(reference_values)
    V  = []; I_br = []; S_br = []; S_inj = [];
    if isfield(reference_values, 'V') 
        V = reference_values.V;    
    end
    
    if isfield(reference_values, 'I') 
        I_br = reference_values.I;
    end
   
    if isfield(reference_values, 'S_br') 
        S_br = reference_values.S_br;
    end
    
    if isfield(reference_values, 'S_inj') 
        S_inj = reference_values.S_inj;
    end
end
    
function [Vmag_nod_idx, Imag_br_idx, Isync_magphase_br_idx, Vsync_magphase_nod_idx, PQ_br_idx, PQ_inj_idx] = extract_meas_indices(meas_indices)
    Vmag_nod_idx = []; Imag_br_idx = []; Isync_magphase_br_idx = []; Vsync_magphase_nod_idx = []; PQ_br_idx = []; PQ_inj_idx = [];
    if isfield(meas_indices, 'Vmag_nod_idx') 
        Vmag_nod_idx = meas_indices.Vmag_nod_idx;
    end

    if isfield(meas_indices, 'Imag_br_idx') 
        Imag_br_idx = meas_indices.Imag_br_idx;
    end    
    
    if isfield(meas_indices, 'Isync_magphase_br_idx') 
        Isync_magphase_br_idx = meas_indices.Isync_magphase_br_idx;
    end

    if isfield(meas_indices, 'Vsync_magphase_nod_idx') 
        Vsync_magphase_nod_idx = meas_indices.Vsync_magphase_nod_idx;
    end
    
    if isfield(meas_indices, 'PQ_br_idx') 
        PQ_br_idx = meas_indices.PQ_br_idx;
    end

    if isfield(meas_indices, 'PQ_inj_idx') 
        PQ_inj_idx = meas_indices.PQ_inj_idx;
    end
end

