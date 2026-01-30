function[V, I, num_iter] = BC_DSSE_rect_trad_3Ph_2polar(num_nodes,nominal_linedata, ...
	node_3ph, node_3ph_reverse, V_abs, V_theta, Ir, Ix, zdata, Znod_state, Gz)


branch_end = nominal_linedata(:,2);
num_branches = length(branch_end);

% 3 Ph topology data
branch_first_A = nominal_linedata(:,3);
branch_first_B = nominal_linedata(:,4);
branch_first_C = nominal_linedata(:,5);
branch_end_A = nominal_linedata(:,6);
branch_end_B = nominal_linedata(:,7);
branch_end_C = nominal_linedata(:,8);
branch_cod_A = nominal_linedata(:,10);
branch_cod_B = nominal_linedata(:,11);
branch_cod_C = nominal_linedata(:,12);

% indici dei rami della sottorete A (B o C) nell'indicizzazione dei rami della rete
A_branches_indices = find(branch_first_A > 0);
B_branches_indices = find(branch_first_B > 0);
C_branches_indices = find(branch_first_C > 0);

num_A_nodes = max(node_3ph(:,1));
num_B_nodes = max(node_3ph(:,2));
num_C_nodes = max(node_3ph(:,3));

% indici dei nodi della whole network that exists in the 3 subnets (at
% indices i, there is the i-th node of A translated in whole net indexing.
A_nodes_indices = node_3ph_reverse(1:num_A_nodes,1);
B_nodes_indices = node_3ph_reverse(1:num_B_nodes,2);
C_nodes_indices = node_3ph_reverse(1:num_C_nodes,3);

% indice dei nodi first di ogni ramo della sottorete, espressi
% con l'indicizzazione della sottorete.
A_branch_first = branch_first_A(A_branches_indices);
B_branch_first = branch_first_B(B_branches_indices);
C_branch_first = branch_first_C(C_branches_indices);
% indice dei nodi end di ogni ramo della sottorete, espressi
% con l'indicizzazione della sottorete.
A_branch_end = branch_end_A(A_branches_indices);
B_branch_end = branch_end_B(B_branches_indices);
C_branch_end = branch_end_C(C_branches_indices);
% numero di rami in ogni sottorete
num_branchesA = sum(branch_cod_A > 0);
num_branchesB = sum(branch_cod_B > 0);
num_branchesC = sum(branch_cod_C > 0);
% numero di variabili di stato = 2 x (somma rami/correnti delle sottoreti)
% + 3 moduli delle tensioni (fasi del primo nodo fissate)
num_state_var = 2*(num_branchesA + num_branchesB + num_branchesC) + 3;     

% correnti dei rami Re e Im delle sottoreti
IrA = Ir(A_branches_indices);
IrB = Ir(B_branches_indices + num_branches);
IrC = Ir(C_branches_indices + 2*(num_branches));
IxA = Ix(A_branches_indices);
IxB = Ix(B_branches_indices + num_branches);
IxC = Ix(C_branches_indices + 2*(num_branches));

type = zdata(:,1);
z = zdata(:,2:4);

% indici dei rami della rete per ogni misura
br = abs(zdata(:,5));   br_sign = sign(zdata(:,5));     br_sign(br_sign==0) = 1;
for i=1:3
    z(:,i) = z(:,i).*br_sign;
end
% indici dei nodi from per ogni misura
fbus = zdata(:,6);
% indici nel vettore delle misure delle misure che riguardano le singole fasi
zAindices = find(node_3ph(fbus,1) ~= 0);
zBindices = find(node_3ph(fbus,2) ~= 0);
zCindices = find(node_3ph(fbus,3) ~= 0);

br1A = br; br1B = br; br1C = br;
br_id = br==0;
brA = find(branch_cod_A,1,'first');
brB = find(branch_cod_B,1,'first');
brC = find(branch_cod_C,1,'first');
br1A(br_id) = brA;
br1B(br_id) = brB;
br1C(br_id) = brC;
zAind = find(branch_cod_A(br1A) == 0);   %%% questa parte é fatta per escludere dall'elenco delle misure quelle relative ai rami inesistenti perché non esiste la corrispondente fase
zBind = find(branch_cod_B(br1B) == 0);
zCind = find(branch_cod_C(br1C) == 0);
fl_measA = size(zAind,1);
fl_measB = size(zBind,1);
fl_measC = size(zCind,1);
if fl_measA > 0
    for i = 1: fl_measA
        zAindices(zAindices==zAind(fl_measA-i+1)) = [];
    end
end
if fl_measB > 0
    for i = 1: fl_measB
        zBindices(zBindices==zBind(fl_measB-i+1)) = [];
    end
end
if fl_measC > 0
    for i = 1: fl_measC
        zCindices(zCindices==zCind(fl_measC-i+1)) = [];
    end
end

% Misure per le varie fasi
zA=z(zAindices,1);
zB=z(zBindices,2);
zC=z(zCindices,3);
% tbus = zdata(:,7);
dev_std = zdata(:,8:10);
sigma2 = dev_std.^2;
sigma2(sigma2<10^-12) =10^-12;
% varianze delle misure per le singole fasi
sigma2A = sigma2(zAindices,1);
sigma2B = sigma2(zBindices,2);
sigma2C = sigma2(zCindices,3);

% vettori di flag corrispondenti alle posizioni nel vettore delle misure per:
vii = (type == 1);      % misure di modulo della tensione
pii = (type == 2);      % misure di potenza attiva iniettata
qii = (type == 3);      % misure di potenza reattiva iniettata
pfi = (type == 4);      % misure di flusso di potenza attiva
qfi = (type == 5);      % misure di flusso di potenza reattiva
ifi = (type == 6);      % misure di modulo della corrente di ramo

% vettori di indici nel vettore delle misure per:
vi = find(vii);         % Index of voltage magnitude measurements..
pidx = find(pii);         % Index of real power injection measurements..
qi = find(qii);         % Index of reactive power injection measurements..
pf = find(pfi);         % Index of real powerflow measurements..
qf = find(qfi);         % Index of reactive powerflow measurements..
iampf = find(ifi);      % Index of current amplitude measurements..

busvi = fbus(vi);       % Node indices (whole network) of Voltage amplitude measurements..
busppi = fbus(pidx);      % Node indices (whole network) of Real Power Injection measurements..
buspqi = fbus(qi);      % Node indices (whole network) of Reactive Power Injection measurements..
buspfi = fbus(pf);
brpf = br(pf);          % Branch indices (whole network) of Active Power flow measurements..
brqf = br(qf);          % Branch indices (whole network) of Reactive Power flow measurements..
briampf = br(iampf);    % Branch indices (whole network) of Current Amplitude flow measurements..

% Impedance matrix to find voltage phasors 2:num_nodes-1 from current phasors (whole network indexing)
% num_nodes-1 x num_branches (che è num_nodes-1)
Rnod = real(Znod_state);
Xnod = imag(Znod_state);

Vr = V_abs.*cos(V_theta);                                                  % calcolo componente reale e immaginaria del vettore delle tensioni nodali 
Vx = V_abs.*sin(V_theta);

% Tensioni dei nodi di ogni singola sottorete, come Re e Im, come Abs e
% Theta e come V complex
VrA = Vr(A_nodes_indices);   VrB = Vr(B_nodes_indices + num_nodes);   VrC = Vr(C_nodes_indices + 2*num_nodes);
VxA = Vx(A_nodes_indices);   VxB = Vx(B_nodes_indices + num_nodes);   VxC = Vx(C_nodes_indices + 2*num_nodes);
V_absA = V_abs(A_nodes_indices);   V_absB = V_abs(B_nodes_indices + num_nodes);   V_absC = V_abs(C_nodes_indices + 2*num_nodes);
V_thetaA = V_theta(A_nodes_indices);   V_thetaB = V_theta(B_nodes_indices + num_nodes);   V_thetaC = V_theta(C_nodes_indices + 2*num_nodes);
VA = complex(VrA,VxA);   VB = complex(VrB,VxB);   VC = complex(VrC,VxC);

% indici dei pesi del WLS, inizializzati per ogni fase a 1:lunghezza del
% vettore di varianza
hidxA = 1 : length(sigma2A); hidxB = 1 : length(sigma2B); hidxC = 1 : length(sigma2C);
vidxA = 1 : length(sigma2A); vidxB = 1 : length(sigma2B); vidxC = 1 : length(sigma2C);
covariancesA = sigma2A; covariancesB = sigma2B; covariancesC = sigma2C;

% indici dei nodi from di sottorete, dei rami sempre di sottorete, numero misure per tipo e sottorete (indicizzazione per ogni sottorete) per le misure di:
% V 
busviA = node_3ph(busvi,1); busviA = busviA(busviA > 0);nviA = length(busviA);
busviB = node_3ph(busvi,2); busviB = busviB(busviB > 0);nviB = length(busviB);
busviC = node_3ph(busvi,3); busviC = busviC(busviC > 0);nviC = length(busviC);
% Pinj
busppiA = node_3ph(busppi,1); busppiA = busppiA(busppiA > 0); npiA = length(busppiA);
busppiB = node_3ph(busppi,2); busppiB = busppiB(busppiB > 0); npiB = length(busppiB);
busppiC = node_3ph(busppi,3); busppiC = busppiC(busppiC > 0); npiC = length(busppiC);
% Qinj
buspqiA = node_3ph(buspqi,1); buspqiA = buspqiA(buspqiA > 0); nqiA = length(buspqiA);
buspqiB = node_3ph(buspqi,2); buspqiB = buspqiB(buspqiB > 0); nqiB = length(buspqiB);
buspqiC = node_3ph(buspqi,3); buspqiC = buspqiC(buspqiC > 0); nqiC = length(buspqiC);
% Pflow
fbuspfA = node_3ph(buspfi,1);
fbuspfB = node_3ph(buspfi,2);
fbuspfC = node_3ph(buspfi,3);
brpfA = branch_cod_A(brpf);     idx_brpfA = find(brpfA==0);
brpfA = brpfA(brpfA > 0);       fbuspfA(idx_brpfA) = [];
npfA = length(brpfA);
brpfB = branch_cod_B(brpf);     idx_brpfB = find(brpfB==0);
brpfB = brpfB(brpfB > 0);       fbuspfB(idx_brpfB) = [];
npfB = length(brpfB);
brpfC = branch_cod_C(brpf);     idx_brpfC = find(brpfC==0);
brpfC = brpfC(brpfC > 0);       fbuspfC(idx_brpfC) = [];
npfC = length(brpfC);
% Qflow
brqfA = branch_cod_A(brqf); brqfA = brqfA(brqfA > 0); nqfA = length(brqfA);
brqfB = branch_cod_B(brqf); brqfB = brqfB(brqfB > 0); nqfB = length(brqfB);
brqfC = branch_cod_C(brqf); brqfC = brqfC(brqfC > 0); nqfC = length(brqfC);
% Iampflow
briampfA = branch_cod_A(briampf); briampfA = briampfA(briampfA > 0); niampfA = length(briampfA);
briampfB = branch_cod_B(briampf); briampfB = briampfB(briampfB > 0); niampfB = length(briampfB);
briampfC = branch_cod_C(briampf); briampfC = briampfC(briampfC > 0); niampfC = length(briampfC);
% Meshes
nmeshA = size(Gz,1)/3; nmeshB = size(Gz,1)/3; nmeshC = size(Gz,1)/3;

% vettori di flag corrispondenti alle posizioni nel vettore delle misure di ogni singola sottorete per:
% Active Power Inj
pii_zdataA = pii(zAindices);  pii_zdataB = pii(zBindices);  pii_zdataC = pii(zCindices);
% Reactive Power Inj
qii_zdataA = qii(zAindices);  qii_zdataB = qii(zBindices);  qii_zdataC = qii(zCindices);
% Active Power Flow Inj
pfi_zdataA = pfi(zAindices);  pfi_zdataB = pfi(zBindices);  pfi_zdataC = pfi(zCindices);
% Reactve Power Flow Inj
qfi_zdataA = qfi(zAindices);  qfi_zdataB = qfi(zBindices);  qfi_zdataC = qfi(zCindices);

%  vettori delle misure di ogni singola sottorete per:
% Active Power Inj
P_injA = zA(pii_zdataA);  P_injB = zB(pii_zdataB);  P_injC = zC(pii_zdataC);
% Reactive Power Inj
Q_injA = zA(qii_zdataA);  Q_injB = zB(qii_zdataB);  Q_injC = zC(qii_zdataC);
% Active Power Flow Inj
P_brA = zA(pfi_zdataA);   P_brB = zB(pfi_zdataB);   P_brC = zC(pfi_zdataC);
% Reactve Power Flow Inj
Q_brA = zA(qfi_zdataA);   Q_brB = zB(qfi_zdataB);   Q_brC = zC(qfi_zdataC);

% measurement start indices per ogni fase (con la sua indicizzazione), in particolare per i seguenti tipi di misure:
% potenza attiva inj
pistartA = nviA + 1;pistartB = nviB + 1;pistartC = nviC + 1;
% potenza reattiva inj
qistartA = pistartA + npiA; qistartB = pistartB + npiB; qistartC = pistartC + npiC; 
% potenza attiva di ramo
pfstartA = qistartA + nqiA; pfstartB = qistartB + nqiB; pfstartC = qistartC + nqiC;
% potenza reattiva di ramo
qfstartA = pfstartA + npfA; qfstartB = pfstartB + npfB; qfstartC = pfstartC + npfC;
% ampiezza corrente di ramo
iampfstartA = qfstartA + nqfA; iampfstartB = qfstartB + nqfB; iampfstartC = qfstartC + nqfC;
% maglie
meshstartA = iampfstartA + niampfA; meshstartB = iampfstartB + niampfB; meshstartC = iampfstartC + niampfC;

zA = [zA; zeros(nmeshA,1)]; zB = [zB; zeros(nmeshB,1)]; zC = [zC; zeros(nmeshC,1)]; 
GzA = Gz(1:nmeshA, :); GzB = Gz(nmeshA+1:nmeshA+nmeshB, :); GzC = Gz(nmeshA+nmeshB+1:end, :);

% Indici di fine dello stato per le fasi
offsetB = num_branchesA + 3;
offsetC = offsetB + num_branchesB;
offsetAx = offsetC + num_branchesC;
offsetBx = offsetAx + num_branchesA;
offsetCx = offsetBx + num_branchesB;

%%% le covarianze calcolate qui di seguito tengono conto del fatto che nel
%%% calcolo delle misure equivalenti di corrente per la fase B e C non si
%%% puó fare la semplificazione che Vr ~= 1 e Vx ~= 0;

%%% Definizione delle sottomatrici Jacobiane costanti

% 2-3) Power Injection Measurements
H2A = zeros(npiA,num_state_var);
H2B = zeros(npiB,num_state_var);
H2C = zeros(npiC,num_state_var);
H3A = zeros(nqiA,num_state_var);
H3B = zeros(nqiB,num_state_var);
H3C = zeros(nqiC,num_state_var);
for i=1:npiA                        % numero misure di potenza attiva iniettata fase A
    m = busppiA(i);                 % indice del nodo misurato di sottorete A (indicizzazione di sottorete)
    toi = (A_branch_end == m);      % vettore dei flag dei rami (di sottorete A) che arrivano al nodo misurato (indicizzazione di sottorete)
    fromi = (A_branch_first == m);  % vettore dei flag dei rami (di sottorete A) che partono dal nodo misurato (indicizzazione di sottorete)
    to = find(toi);                 % vettore degli indici dei rami (di sottorete A) che arrivano al nodo misurato (indicizzazione di sottorete)
    from = find(fromi);             % vettore degli indici dei rami (di sottorete A) che partono dal nodo misurato (indicizzazione di sottorete)
    H2A(i, to+3) = 1;               % + 1 per i rami che arrivano (erogano) 
    H2A(i, from+3) = -1;            % - 1 per i rami che partono (iniettano)
    H3A(i, to + offsetAx) = 1;      % + 1 per i rami che arrivano (erogano)
    H3A(i, from + offsetAx) = -1;   % - 1 per i rami che partono (iniettano)
end
for i=1:npiB                        % numero misure di potenza attiva iniettata fase B
    m = busppiB(i);                 % indice del nodo misurato di sottorete B (indicizzazione di sottorete)
    toi = (B_branch_end == m);      % vettore dei flag dei rami (di sottorete B) che arrivano al nodo misurato (indicizzazione di sottorete)
    fromi = (B_branch_first == m);  % vettore dei flag dei rami (di sottorete B) che partono dal nodo misurato (indicizzazione di sottorete)
    to = find(toi);                 % vettore degli indici dei rami (di sottorete B) che arrivano al nodo misurato (indicizzazione di sottorete)
    from = find(fromi);             % vettore degli indici dei rami (di sottorete B) che partono dal nodo misurato (indicizzazione di sottorete)
    H2B(i, to + offsetB) = 1;       % + 1 per i rami che arrivano (erogano) 
    H2B(i, from + offsetB) = -1;    % - 1 per i rami che partono (iniettano)
    H3B(i, to + offsetBx) = 1;      % + 1 per i rami che arrivano (erogano)
    H3B(i, from + offsetBx) = -1;   % - 1 per i rami che partono (iniettano)
    idx_p = nviB + i;
    idx_q = nviB + npiB + i;
    rot_mat = [cos(2*pi/3) , sin(2*pi/3); sin(2*pi/3), -cos(2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2B(idx_p), 0; 0, sigma2B(idx_q)] * rot_mat';
    covariancesB(idx_p) = covar_pq_eq(1,1);
    covariancesB(idx_q) = covar_pq_eq(2,2);
    hidxB = [hidxB, idx_p, idx_q];
    vidxB = [vidxB, idx_q, idx_p];
    covariancesB = [covariancesB; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end
for i=1:npiC                        % numero misure di potenza attiva iniettata fase C
    m = busppiC(i);                 % indice del nodo misurato di sottorete C (indicizzazione di sottorete)
    toi = (C_branch_end == m);      % vettore dei flag dei rami (di sottorete C) che arrivano al nodo misurato (indicizzazione di sottorete)
    fromi = (C_branch_first == m);  % vettore dei flag dei rami (di sottorete C) che partono dal nodo misurato (indicizzazione di sottorete)
    to = find(toi);                 % vettore degli indici dei rami (di sottorete C) che arrivano al nodo misurato (indicizzazione di sottorete)
    from = find(fromi);             % vettore degli indici dei rami (di sottorete C) che partono dal nodo misurato (indicizzazione di sottorete)
    H2C(i, to + offsetC) = 1;       % + 1 per i rami che arrivano (erogano) 
    H2C(i, from + offsetC) = -1;    % - 1 per i rami che partono (iniettano)
    H3C(i, to + offsetCx) = 1;      % + 1 per i rami che arrivano (erogano)
    H3C(i, from + offsetCx) = -1;   % - 1 per i rami che partono (iniettano)
    idx_p = nviC + i;
    idx_q = nviC + npiC + i;
    rot_mat = [cos(-2*pi/3) , sin(-2*pi/3); sin(-2*pi/3), -cos(-2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2C(idx_p), 0; 0, sigma2C(idx_q)] * rot_mat';
    covariancesC(idx_p) = covar_pq_eq(1,1);
    covariancesC(idx_q) = covar_pq_eq(2,2);
    hidxC = [hidxC, idx_p, idx_q];
    vidxC = [vidxC, idx_q, idx_p];
    covariancesC = [covariancesC; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end

% 4-5) Power Flow Measurements
H4A = zeros(npfA,num_state_var);
H4B = zeros(npfB,num_state_var);
H4C = zeros(npfC,num_state_var);
H5A = zeros(nqfA,num_state_var);
H5B = zeros(nqfB,num_state_var);
H5C = zeros(nqfC,num_state_var);
for i=1:npfA                         % numero misure di potenza attiva iniettata fase A
    idx = brpfA(i);                  % indice del ramo misurato di sottorete A (indicizzazione di sottorete)
    H4A(i, idx + 3) = 1;             % reale
    H5A(i, idx + offsetAx) = 1;      % imag
end
for i=1:npfB                         % numero misure di potenza attiva iniettata fase B
    idx = brpfB(i);                  % indice del ramo misurato di sottorete B (indicizzazione di sottorete)
    H4B(i, idx + offsetB) = 1;       % reale
    H5B(i, idx + offsetBx) = 1;      % imag
    idx_p = pfstartB-1 + i;
    idx_q = qfstartB-1 + i;
    rot_mat = [cos(2*pi/3) , sin(2*pi/3); sin(2*pi/3), -cos(2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2B(idx_p), 0; 0, sigma2B(idx_q)] * rot_mat';
    covariancesB(idx_p) = covar_pq_eq(1,1);
    covariancesB(idx_q) = covar_pq_eq(2,2);
    hidxB = [hidxB, idx_p, idx_q];
    vidxB = [vidxB, idx_q, idx_p];
    covariancesB = [covariancesB; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end
for i=1:npfC                         % numero misure di potenza attiva iniettata fase C
    idx = brpfC(i);                  % indice del ramo misurato di sottorete C (indicizzazione di sottorete)
    H4C(i, idx + offsetC) = 1;       % reale
    H5C(i, idx + offsetCx) = 1;      % imag
    idx_p = pfstartC-1 + i;
    idx_q = qfstartC-1 + i;
    rot_mat = [cos(-2*pi/3) , sin(-2*pi/3); sin(-2*pi/3), -cos(-2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2C(idx_p), 0; 0, sigma2C(idx_q)] * rot_mat';
    covariancesC(idx_p) = covar_pq_eq(1,1);
    covariancesC(idx_q) = covar_pq_eq(2,2);
    hidxC = [hidxC, idx_p, idx_q];
    vidxC = [vidxC, idx_q, idx_p];
    covariancesC = [covariancesC; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end

% 11) Mesh Virtual Measurements
H11A = zeros(nmeshA,num_state_var);
H11B = zeros(nmeshB,num_state_var);
H11C = zeros(nmeshC,num_state_var);
if size(GzA,1) > 0
    for i=1:nmeshA
        H11A(i, :) = [zeros(1,3),GzA(i,:)];
        idx = meshstartA + i - 1;
        hidxA = [hidxA, idx];
        vidxA = [vidxA, idx];
        covariancesA = [covariancesA; 10^-12];
    end
else 
    H11A = [];
end
if size(GzB,1) > 0
    for i=1:nmeshB
        H11B(i, :) = [zeros(1,3),GzB(i,:)];
        idx = meshstartB + i - 1;
        hidxB = [hidxB, idx];
        vidxB = [vidxB, idx];
        covariancesB = [covariancesB; 10^-12];
    end
else
    H11B = [];
end
if size(GzC,1) > 0
    for i=1:nmeshC
        H11C(i, :) = [zeros(1,3),GzC(i,:)];
        idx = meshstartC + i - 1;
        hidxC = [hidxC, idx];
        vidxC = [vidxC, idx];
        covariancesC = [covariancesC; 10^-12];
    end
else
    H11C = [];
end

% Matrice Pesi
% Phase A
WA = sparse(hidxA, vidxA, covariancesA);
% Phase B
WB = sparse(hidxB, vidxB, covariancesB);
% Phase C
WC = sparse(hidxC, vidxC, covariancesC);
% Tutte le 3 fasi 
W = blkdiag(WA, WB, WC);
W = full(W);
% Inversa della W
iW = W\eye(size(W));
iW = sparse(iW);

% Initial state
State = [V_absA(1);V_absB(1);V_absC(1);IrA;IrB;IrC;IxA;IxB;IxC];

num_iter = 0;
epsilonV = 5;
while (epsilonV > 10^-7 && num_iter <= 100)
    
    % Save State and Calculate Equivalent Currents from prev estimation (for residual calculation)
    % Phase A
    V_abs_last_stepA = V_absA;
    V_theta_last_stepA = V_thetaA;
    StateV_lastA = [V_theta_last_stepA; V_abs_last_stepA];                 % stato equivalente a quello degli stimatori di tensione
    
    % Phase B
    V_abs_last_stepB = V_absB;
    V_theta_last_stepB = V_thetaB;
    StateV_lastB = [V_theta_last_stepB; V_abs_last_stepB];                 % stato equivalente a quello degli stimatori di tensione
    
    % Phase C
    V_abs_last_stepC = V_absC;
    V_theta_last_stepC = V_thetaC;
    StateV_lastC = [V_theta_last_stepC; V_abs_last_stepC];                 % stato equivalente a quello degli stimatori di tensione
    
    % Phase A
    % Correnti equivalenti di Ramo 
    Ir_brA = (P_brA.*VrA(fbuspfA)+Q_brA.*VxA(fbuspfA))./(V_absA(fbuspfA).^2);
    Ix_brA = (P_brA.*VxA(fbuspfA)-Q_brA.*VrA(fbuspfA))./(V_absA(fbuspfA).^2);
    zA(pfstartA : pfstartA + npfA -1) = Ir_brA;
    zA(qfstartA : qfstartA + nqfA -1) = Ix_brA;
    % Correnti equivalenti erogate di Nodo (dalle Potenze iniettate) (tutti i nodi da 2 a num nodi fase C
    Ir_injA = (P_injA.*VrA(2:num_A_nodes)+Q_injA.*VxA(2:num_A_nodes))./(V_absA(2:num_A_nodes).^2);
    Ix_injA = (P_injA.*VxA(2:num_A_nodes)-Q_injA.*VrA(2:num_A_nodes))./(V_absA(2:num_A_nodes).^2);
    zA(pistartA : pistartA + npiA -1) = -Ir_injA;          % negativo perchè per convenzione le correnti erogate verso i carichi sono negative e quindi 
    zA(qistartA : qistartA + nqiA -1) = -Ix_injA;          % il segno meno mi serve per bilanciare la relazione relativa ai residui
    
    % Phase B
    % Correnti equivalenti di Ramo 
    Ir_brB = (P_brB.*VrB(fbuspfB)+Q_brB.*VxB(fbuspfB))./(V_absB(fbuspfB).^2);
    Ix_brB = (P_brB.*VxB(fbuspfB)-Q_brB.*VrB(fbuspfB))./(V_absB(fbuspfB).^2);
    zB(pfstartB : pfstartB + npfB -1) = Ir_brB;
    zB(qfstartB : qfstartB + nqfB -1) = Ix_brB;
    % Correnti equivalenti erogate di Nodo (dalle Potenze iniettate) (tutti i nodi da 2 a num nodi fase C
    Ir_injB = (P_injB.*VrB(2:num_B_nodes)+Q_injB.*VxB(2:num_B_nodes))./(V_absB(2:num_B_nodes).^2);
    Ix_injB = (P_injB.*VxB(2:num_B_nodes)-Q_injB.*VrB(2:num_B_nodes))./(V_absB(2:num_B_nodes).^2);
    zB(pistartB : pistartB + npiB -1) = -Ir_injB;          % negativo perchè per convenzione le correnti erogate verso i carichi sono negative e quindi 
    zB(qistartB : qistartB + nqiB -1) = -Ix_injB;          % il segno meno mi serve per bilanciare la relazione relativa ai residui
    
    % Phase C
    % Correnti equivalenti di Ramo 
    Ir_brC = (P_brC.*VrC(fbuspfC)+Q_brC.*VxC(fbuspfC))./(V_absC(fbuspfC).^2);
    Ix_brC = (P_brC.*VxC(fbuspfC)-Q_brC.*VrC(fbuspfC))./(V_absC(fbuspfC).^2);
    zC(pfstartC : pfstartC + npfC -1) = Ir_brC;
    zC(qfstartC : qfstartC + nqfC -1) = Ix_brC;
    % Correnti equivalenti erogate di Nodo (dalle Potenze iniettate) (tutti i nodi da 2 a num nodi fase C)
    Ir_injC = (P_injC.*VrC(2:num_C_nodes)+Q_injC.*VxC(2:num_C_nodes))./(V_absC(2:num_C_nodes).^2);
    Ix_injC = (P_injC.*VxC(2:num_C_nodes)-Q_injC.*VrC(2:num_C_nodes))./(V_absC(2:num_C_nodes).^2);
    zC(pistartC : pistartC + npiC -1) = -Ir_injC;          % negativo perchè per convenzione le correnti erogate verso i carichi sono negative e quindi 
    zC(qistartC : qistartC + nqiC -1) = -Ix_injC;          % il segno meno mi serve per bilanciare la relazione relativa ai residui
    
    
    %%% 1) Voltage Magnitude Measurements   
    %Phase A
    h1A = V_absA(busviA,1); 
    H1A = zeros(nviA,num_state_var);
    for i=1:nviA                  % numero misure di modulo V fase A
        index = busviA(i);        % indice del nodo di sottorete A (indicizzazione di sottorete)
        idx = index-1;
        H1A(i, 1) = cos(V_thetaA(1) - V_thetaA(index));
        if index == 1             
           continue
        else
            H1A(i, 4:offsetAx) = - Rnod(idx,:)*cos(V_thetaA(index)) - Xnod(idx,:)*sin(V_thetaA(index));
            H1A(i, offsetAx+1:end) = - Rnod(idx,:)*sin(V_thetaA(index)) + Xnod(idx,:)*cos(V_thetaA(index));
        end
    end
    %Phase B
    h1B = V_absB(busviB,1); 
    H1B = zeros(nviB,num_state_var);
    for i=1:nviB                 % numero misure di modulo V fase B
        index = busviB(i);       % indice del nodo di sottorete B (indicizzazione di sottorete)
        idx = index + num_A_nodes - 2;
        H1B(i, 2) = cos(V_thetaB(1) - V_thetaB(index));
        if index == 1            % se sto misurando il nodo 1 fase B che è nello stato
            continue
        else
            H1B(i, 4:offsetAx) = - Rnod(idx,:)*cos(V_thetaB(index)) - Xnod(idx,:)*sin(V_thetaB(index));
            H1B(i, offsetAx+1:end) = - Rnod(idx,:)*sin(V_thetaB(index)) + Xnod(idx,:)*cos(V_thetaB(index));
        end
    end
    %Phase C
    h1C = V_absC(busviC,1);
    H1C = zeros(nviC,num_state_var);
    for i=1:nviC    % numero misure di modulo V fase C
        index = busviC(i);        % indice del nodo di sottorete C (indicizzazione di sottorete)
        idx = index + num_A_nodes + num_B_nodes - 3;
        H1C(i, 3) = cos(V_thetaC(1) - V_thetaC(index));
        if index == 1             % se sto misurando il nodo 1 fase C che è nello stato
            continue
        else
            H1C(i, 4:offsetAx) = - Rnod(idx,:)*cos(V_thetaC(index)) - Xnod(idx,:)*sin(V_thetaC(index));
            H1C(i, offsetAx+1:end) = - Rnod(idx,:)*sin(V_thetaC(index)) + Xnod(idx,:)*cos(V_thetaC(index));
        end
    end
    
    %%% 2-3) Power Injection Measurements   (OK, come indicato nel paper di Baran Kelley)
    h2A = H2A*State;
    h2B = H2B*State;
    h2C = H2C*State;
    h3A = H3A*State;
    h3B = H3B*State;
    h3C = H3C*State;   
    
    %%% 4-5) Power Flow Measurements   (OK, come indicato nel paper di Baran Kelley)
    h4A = H4A*State;
    h4B = H4B*State;
    h4C = H4C*State;
    h5A = H5A*State;
    h5B = H5B*State;
    h5C = H5C*State;
    
    %%% 6) Current Magnitude Measurements   (OK, come indicato nel paper di Baran Kelley)
    h6A = sqrt(IrA(briampfA).^2+IxA(briampfA).^2);
    H6A = zeros(niampfA,num_state_var);
    I_Total_Complex = complex(IrA(briampfA),IxA(briampfA));
    phi_Imag_br = angle(I_Total_Complex);
    for i=1:niampfA                         % numero misure di potenza attiva iniettata fase A
        idx = briampfA(i);              % indice del ramo misurato di sottorete A (indicizzazione di sottorete)                    
        H6A(i, idx + 3) = cos(phi_Imag_br(i));
        H6A(i, idx + offsetAx) = sin(phi_Imag_br(i));
    end
    h6B = sqrt(IrB(briampfB).^2+IxB(briampfB).^2);
    H6B = zeros(niampfB,num_state_var);
    I_Total_Complex = complex(IrB(briampfB),IxB(briampfB));
    phi_Imag_br = angle(I_Total_Complex);
    for i=1:niampfB                         % numero misure di potenza attiva iniettata fase B
        idx = briampfB(i);              % indice del ramo misurato di sottorete B (indicizzazione di sottorete)
        H6B(i, idx + offsetB) = cos(phi_Imag_br(i));
        H6B(i, idx + offsetBx) = sin(phi_Imag_br(i));
    end
    h6C = sqrt(IrC(briampfC).^2+IxC(briampfC).^2);
    H6C = zeros(niampfC,num_state_var);
    I_Total_Complex = complex(IrC(briampfC),IxC(briampfC));
    phi_Imag_br = angle(I_Total_Complex);
    for i=1:niampfC                         % numero misure di potenza attiva iniettata fase C
        idx = briampfC(i);              % indice del ramo misurato di sottorete C (indicizzazione di sottorete)
        H6C(i, idx + offsetC) = cos(phi_Imag_br(i));
        H6C(i, idx + offsetCx) = sin(phi_Imag_br(i));
    end
    
    %%% 11) Mesh Virtual Measurements
    if size(GzA,1) > 0
        h11A = H11A*State;
    else
        h11A = [];
    end
    if size(GzB,1) > 0
        h11B = H11B*State;
    else
        h11B = [];
    end
    if size(GzC,1) > 0
        h11C = H11C*State;
    else
        h11C = [];
    end
        
    %%%  fine misure
    
    % Composizione vettore misure e stima
    % Phase A
    hA = [h1A; h2A; h3A; h4A; h5A; h6A; h11A];
    HA = [H1A; H2A; H3A; H4A; H5A; H6A; H11A];
    % Phase B
    hB = [h1B; h2B; h3B; h4B; h5B; h6B; h11B];
    HB = [H1B; H2B; H3B; H4B; H5B; H6B; H11B];
    % Phase C
    hC = [h1C; h2C; h3C; h4C; h5C; h6C; h11C];
    HC = [H1C; H2C; H3C; H4C; H5C; H6C; H11C];
    
    z = [zA; zB; zC];
    h = [hA; hB; hC];
    r = z-h;
    H = [HA; HB; HC];   % Jacobiano completo
    H = sparse(H);
    g = H'*iW*r;
    Gm = H'*iW*H;

    Delta_X = Gm\g;                       % WLS di stima dello stato
    State = State + Delta_X;              % aggiornamento variabili di stato
    
    % Phase A
    VA(1) = State(1) * exp(1i*V_thetaA(1));
    IrA = State(4 : offsetB);
    IxA = State(offsetAx + 1 : offsetBx);
    % Phase B
    VB(1) = State(2)*exp(1i * V_thetaB(1));
    IrB = State(offsetB + 1 : offsetC);
    IxB = State(offsetBx + 1 : offsetCx);
    % Phase C
    VC(1) = State(3) * exp(1i * V_thetaC(1));
    IrC = State(offsetC + 1 : offsetAx);
    IxC = State(offsetCx + 1 : end);
    
    % FORWARD SWEEP tri-fase ATTENZIONE PER RADIALI
    IcomplexA = complex(IrA, IxA);
    IcomplexB = complex(IrB, IxB);
    IcomplexC = complex(IrC, IxC);
    I_ABC_Complex = [IcomplexA; IcomplexB; IcomplexC];
    % 
    % VA(2:end) = VA(1) - Znod_state(1:num_A_nodes - 1,:) * I_ABC_Complex;
    % VB(2:end) = VB(1) - Znod_state(num_A_nodes : num_A_nodes + num_B_nodes-2, :)*I_ABC_Complex;
    % VC(2:end) = VC(1) - Znod_state(num_A_nodes + num_B_nodes - 1: end, :)*I_ABC_Complex;
    % 
    % %Vincenzo, 
    % V = zeros(size(Vr));    
    % V(A_nodes_indices, 1) = VA;
    % V(B_nodes_indices+num_nodes, 1) = VB;
    % V(C_nodes_indices+2*num_nodes, 1) = VC;

    %% Vincenzo 14/07/2025
    % Riscritte le righe sopra in forma matriciale e verificatta la
    % correttezza V3Ph = V1_3ph - Z_node_state_3ph * Icomplex =
    % Trasformazione * [V1_3ph; I_ABC_Complex]
    % Attenzione, della tensione al nodo 1 stimiamo SOLO i moduli (A, B, C)
    V1_3ph = [VA(1); VB(1); VC(1)];
    num_ABC_nodes = num_A_nodes + num_B_nodes + num_C_nodes;
    maskV1A = zeros(num_A_nodes, 3); maskV1A(:, 1) = 1;
    maskV1B = zeros(num_B_nodes, 3); maskV1B(:, 2) = 1;
    maskV1C = zeros(num_C_nodes, 3); maskV1C(:, 3) = 1;
    MaskV1 = [maskV1A; maskV1B; maskV1C];  
   
    Znod_state_A = Znod_state(1:num_A_nodes - 1,:);
    Znod_state_B = Znod_state(num_A_nodes : num_A_nodes + num_B_nodes - 2,:);
    Znod_state_C = Znod_state(num_A_nodes + num_B_nodes - 1: end,:);
    zeros_row = zeros(1, num_ABC_nodes - 3);

    TrSC2V_Complex =  [MaskV1, [zeros_row; -Znod_state_A; ...
                               zeros_row; -Znod_state_B; ...
                               zeros_row; -Znod_state_C]];
    StateComplex = [V1_3ph; I_ABC_Complex];
    V_ABC_Complex = TrSC2V_Complex * StateComplex;

    Vcomplex = zeros(size(Vr));    
    Vcomplex(A_nodes_indices, 1) = V_ABC_Complex(1 : num_A_nodes);
    Vcomplex(B_nodes_indices + num_nodes, 1) = V_ABC_Complex(num_A_nodes + 1 : num_A_nodes + num_B_nodes);
    Vcomplex(C_nodes_indices + 2 * num_nodes, 1) = V_ABC_Complex(num_A_nodes + num_B_nodes + 1 : end);
    
    Vr = real(Vcomplex);
    Vx = imag(Vcomplex);
    VrA = Vr(A_nodes_indices);VrB = Vr(B_nodes_indices + num_nodes);VrC = Vr(C_nodes_indices + 2*num_nodes);
    VxA = Vx(A_nodes_indices);VxB = Vx(B_nodes_indices + num_nodes);VxC = Vx(C_nodes_indices + 2*num_nodes);
    V_abs = abs(Vcomplex);
    V_theta = angle(Vcomplex);
    V_absA = V_abs(A_nodes_indices);V_absB = V_abs(B_nodes_indices + num_nodes);V_absC = V_abs(C_nodes_indices + 2*num_nodes);
    V_thetaA = V_theta(A_nodes_indices);V_thetaB = V_theta(B_nodes_indices + num_nodes);V_thetaC = V_theta(C_nodes_indices + 2*num_nodes);
    
    num_iter = num_iter + 1;
    
    % verifica epsilon di convergenza
    % Phase A
    StateVA = [V_thetaA; V_absA];  % attenzione alle fasi ATTENZIONE
    Delta_VA = StateVA - StateV_lastA;
    % Phase B
    StateVB = [V_thetaB; V_absB];  % attenzione alle fasi ATTENZIONE
    Delta_VB = StateVB - StateV_lastB;
    % Phase C
    StateVC = [V_thetaC; V_absC];  % attenzione alle fasi ATTENZIONE
    Delta_VC = StateVC - StateV_lastC;
    Delta_V = [Delta_VA; Delta_VB; Delta_VC];

    epsilonV = max(abs(Delta_V));
end

Ir(A_branches_indices) = IrA;Ir(B_branches_indices + num_branches) = IrB;Ir(C_branches_indices + 2*num_branches) = IrC;
Ix(A_branches_indices) = IxA;Ix(B_branches_indices + num_branches) = IxB;Ix(C_branches_indices + 2*num_branches) = IxC;

%%Vincenzo 2025 EstimationUncertainty
% 

%Covarianza della stima dello stato 
%State = [V_absA(1);V_absB(1);V_absC(1);IrA;IrB;IrC;IxA;IxB;IxC];
Sigma_State = Gm \ eye(size(Gm));

%StateComplex = [V_absA(1)*exp(1i*V_thetaA(1),
%V_absB(1)*exp(1i*V_thetaB(1)),
%V_absC(1)*exp(1i*V_thetaC(1)),IA_complex,IB_complex, IC_complex
%VABC = TrSC2V_Complex * StateComplex;


Znod_state_A_real = real(Znod_state_A);
Znod_state_A_imag = imag(Znod_state_A);
Znod_state_B_real = real(Znod_state_B);
Znod_state_B_imag = imag(Znod_state_B);
Znod_state_C_real = real(Znod_state_C);
Znod_state_C_imag = imag(Znod_state_C);

% V_rect = [VABC_real; VABC_Imag] = Transformation_State2V * State;
Transformation_State2V = [%Reale part
                          cos(V_thetaA(1)) * maskV1A, [zeros_row, zeros_row; 
                                                      Znod_state_A_real, -Znod_state_A_imag];
                          cos(V_thetaB(1)) * maskV1B, [zeros_row, zeros_row; 
                                                      Znod_state_B_real, -Znod_state_B_imag];
                          cos(V_thetaC(1)) * maskV1C, [zeros_row, zeros_row; 
                                                      Znod_state_C_real, -Znod_state_C_imag];
                          % immaginary part
                          sin(V_thetaA(1)) * maskV1A, [zeros_row, zeros_row; 
                                                      Znod_state_A_imag, Znod_state_A_real];
                          sin(V_thetaB(1)) * maskV1B, [zeros_row, zeros_row; 
                                                      Znod_state_B_imag, Znod_state_B_real];
                          sin(V_thetaC(1)) * maskV1C, [zeros_row, zeros_row; 
                                                      Znod_state_C_imag, Znod_state_C_real]];
%Sigma_V_rect = Transformation_State2V * Sigma_State *
%Transformation_State2V';
Sigma_V_rect = Transformation_State2V * Sigma_State * Transformation_State2V';
%%
%% J = [cos -r sin;
%%        sin  r cos];
%Sigma_RX = J * Sigma_polar * J';
%Sigma_polar = iJ * Sigma_RX * iJ';
VABC_Module = abs(V_ABC_Complex);
VABC_Angle = angle(V_ABC_Complex);

J_V = [diag(cos(VABC_Angle)), -diag(VABC_Module .* sin(VABC_Angle));
    diag(sin(VABC_Angle)),   diag(VABC_Module .* cos(VABC_Angle))];
iJ_V = J_V \ eye(size(J_V));

Sigma_V_polar = iJ_V * Sigma_V_rect * iJ_V';
std_V = real(reshape(sqrt(diag(Sigma_V_polar)), [], 2));

%% Current
Sigma_I_rect = Sigma_State(4 : end, 4 : end);

IABC_Module = abs(I_ABC_Complex);
IABC_Angle = angle(I_ABC_Complex);

J_I = [diag(cos(IABC_Angle)), -diag(IABC_Module .* sin(IABC_Angle));
    diag(sin(IABC_Angle)),   diag(IABC_Module .* cos(IABC_Angle))];
iJ_I = J_I \ eye(size(J_I));

Sigma_I_polar = iJ_I * Sigma_I_rect * iJ_I';
std_I = real(reshape(sqrt(diag(Sigma_I_polar)), [], 2));

coverage_factor = 3;

V_Total_Complex = complex(Vr, Vx);
I_Total_Complex = complex(Ir, Ix);

nodes_indices = [A_nodes_indices; B_nodes_indices + num_nodes; C_nodes_indices + 2 * num_nodes];
branches_indices = [A_branches_indices; B_branches_indices + num_branches; C_branches_indices + 2 * num_branches];

V.Module = nan(size(V_Total_Complex));
V.Angle = nan(size(V_Total_Complex));
I.Module = nan(size(I_Total_Complex));
I.Angle = nan(size(I_Total_Complex));

V.Module(nodes_indices) = abs(V_Total_Complex(nodes_indices));
V.Angle(nodes_indices) = angle(V_Total_Complex(nodes_indices));
I.Module(branches_indices) = abs(I_Total_Complex(branches_indices));
I.Angle(branches_indices) = angle(I_Total_Complex(branches_indices));

V.Unc.Module = nan(size(V_Total_Complex));
V.Unc.Angle  = nan(size(V_Total_Complex));
I.Unc.Module = nan(size(I_Total_Complex));
I.Unc.Angle  = nan(size(I_Total_Complex));

V.Unc.Module(nodes_indices) = std_V(:, 1) * coverage_factor;
V.Unc.Angle(nodes_indices)  = std_V(:, 2) * coverage_factor;
I.Unc.Module(branches_indices) = std_I(:, 1) * coverage_factor;
I.Unc.Angle(branches_indices)  = std_I(:, 2) * coverage_factor;





