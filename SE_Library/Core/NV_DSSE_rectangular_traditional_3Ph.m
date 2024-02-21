function[Vr, Vx, Ir, Ix, num_iter] = NV_DSSE_rectangular_traditional_3Ph(num_nodes,nominal_linedata, ...
	node_3ph, node_3ph_reverse, V_abs, V_theta, I_abs, I_theta, zdata, Y_state, Y_tot)

branch_first = nominal_linedata(:,1);
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

num_A_nodes = max(node_3ph(:,1));
num_B_nodes = max(node_3ph(:,2));
num_C_nodes = max(node_3ph(:,3));

% indici dei nodi della whole network that exists in the 3 subnets (at
% indices i, there is the i-th node of A translated in whole net indexing.
A_nodes_indices = node_3ph_reverse(1:num_A_nodes,1);
B_nodes_indices = node_3ph_reverse(1:num_B_nodes,2);
C_nodes_indices = node_3ph_reverse(1:num_C_nodes,3);

% numero di variabili di stato = 2 x (somma nodi della rete)
num_state_var = 2*(num_A_nodes + num_B_nodes + num_C_nodes);     

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

% Admittance matrix 
G_state = real(Y_state);  % parte reale matrice di ammettenza
B_state = imag(Y_state);  % parte immaginaria matrice di ammettenza
Yabs_state = abs(Y_state);
Yphase_state = angle(Y_state);

Vr = V_abs.*cos(V_theta);                                                  % calcolo componente reale e immaginaria del vettore delle tensioni nodali 
Vx = V_abs.*sin(V_theta);

% Tensioni dei nodi di ogni singola sottorete, come Re e Im, come Abs e
% Theta e come V complex
VrA = Vr(A_nodes_indices);VrB = Vr(B_nodes_indices + num_nodes);VrC = Vr(C_nodes_indices + 2*num_nodes);
Vr_red = [VrA; VrB; VrC];
VxA = Vx(A_nodes_indices);VxB = Vx(B_nodes_indices + num_nodes);VxC = Vx(C_nodes_indices + 2*num_nodes);
Vx_red = [VxA; VxB; VxC];
V_absA = V_abs(A_nodes_indices);V_absB = V_abs(B_nodes_indices + num_nodes);V_absC = V_abs(C_nodes_indices + 2*num_nodes);
V_thetaA = V_theta(A_nodes_indices);V_thetaB = V_theta(B_nodes_indices + num_nodes);V_thetaC = V_theta(C_nodes_indices + 2*num_nodes);

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
tbuspfA = branch_end_A(brpf); tbuspfA = tbuspfA(tbuspfA > 0);       % andrebbe rivista se si usassero le potenze in direzione opposta a quella del ramo
tbuspfB = branch_end_B(brpf); tbuspfB = tbuspfB(tbuspfB > 0);
tbuspfC = branch_end_C(brpf); tbuspfC = tbuspfC(tbuspfC > 0);
% Qflow
brqfA = branch_cod_A(brqf); brqfA = brqfA(brqfA > 0); nqfA = length(brqfA);
brqfB = branch_cod_B(brqf); brqfB = brqfB(brqfB > 0); nqfB = length(brqfB);
brqfC = branch_cod_C(brqf); brqfC = brqfC(brqfC > 0); nqfC = length(brqfC);
% Iampflow
fbusiampfA = branch_first_A(briampf); fbusiampfA = fbusiampfA(fbusiampfA > 0);
fbusiampfB = branch_first_B(briampf); fbusiampfB = fbusiampfB(fbusiampfB > 0);
fbusiampfC = branch_first_C(briampf); fbusiampfC = fbusiampfC(fbusiampfC > 0);
tbusiampfA = branch_end_A(briampf); tbusiampfA = tbusiampfA(tbusiampfA > 0);
tbusiampfB = branch_end_B(briampf); tbusiampfB = tbusiampfB(tbusiampfB > 0);
tbusiampfC = branch_end_C(briampf); tbusiampfC = tbusiampfC(tbusiampfC > 0);
briampfA = branch_cod_A(briampf); briampfA = briampfA(briampfA > 0); niampfA = length(briampfA);
briampfB = branch_cod_B(briampf); briampfB = briampfB(briampfB > 0); niampfB = length(briampfB);
briampfC = branch_cod_C(briampf); briampfC = briampfC(briampfC > 0); niampfC = length(briampfC);

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
% iampfstartA = qfstartA + nqfA; iampfstartB = qfstartB + nqfB; iampfstartC = qfstartC + nqfC;

% Indice di fine dello stato della fase A (Inizio B - 1)
offsetB = num_A_nodes;                      
offsetC = num_A_nodes + num_B_nodes;
offsetAx = offsetC + num_C_nodes;
offsetBx = offsetAx + num_A_nodes -1;
offsetCx = offsetBx + num_B_nodes -1;
no_col = [offsetCx+3,offsetBx+2,offsetAx+1];

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
for i=1:npiA
    mA = busppiA(i);
    H2A(i, :) = [-G_state(mA,:),B_state(mA,:)];
    H3A(i, :) = [-B_state(mA,:),-G_state(mA,:)];
end
for i=1:npiB
    mB = busppiB(i) + offsetB;
    H2B(i, :) = [-G_state(mB,:),B_state(mB,:)];  
    H3B(i, :) = [-B_state(mB,:),-G_state(mB,:)]; 
    idx_p = nviB + i;
    idx_q = nviB + npiB +i;
    rot_mat = [cos(2*pi/3) , sin(2*pi/3); sin(2*pi/3), -cos(2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2B(idx_p), 0; 0, sigma2B(idx_q)] * rot_mat';
    covariancesB(idx_p) = covar_pq_eq(1,1);
    covariancesB(idx_q) = covar_pq_eq(2,2);
    hidxB = [hidxB, idx_p, idx_q];
    vidxB = [vidxB, idx_q, idx_p];
    covariancesB = [covariancesB; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end
for i=1:npiC
    mC = busppiC(i) + offsetC;
    H2C(i, :) = [-G_state(mC,:),B_state(mC,:)];  
    H3C(i, :) = [-B_state(mC,:),-G_state(mC,:)];
    idx_p = nviC + i;
    idx_q = nviC + npiC +i;
    rot_mat = [cos(4*pi/3) , sin(4*pi/3); sin(4*pi/3), -cos(4*pi/3)];
    covar_pq_eq = rot_mat * [sigma2C(idx_p), 0; 0, sigma2C(idx_q)] * rot_mat';
    covariancesC(idx_p) = covar_pq_eq(1,1);
    covariancesC(idx_q) = covar_pq_eq(2,2);
    hidxC = [hidxC, idx_p, idx_q];
    vidxC = [vidxC, idx_q, idx_p];
    covariancesC = [covariancesC; covar_pq_eq(1,2); covar_pq_eq(2,1)];
end
H2A(:,offsetB+1) = H2A(:,offsetB+1) + H2A(:,offsetBx+2)*tan(2*pi/3);
H2A(:,offsetC+1) = H2A(:,offsetC+1) + H2A(:,offsetCx+3)*tan(4*pi/3);
H2A(:,no_col) = [];
H2B(:,offsetB+1) = H2B(:,offsetB+1) + H2B(:,offsetBx+2)*tan(2*pi/3);
H2B(:,offsetC+1) = H2B(:,offsetC+1) + H2B(:,offsetCx+3)*tan(4*pi/3);
H2B(:,no_col) = [];
H2C(:,offsetB+1) = H2C(:,offsetB+1) + H2C(:,offsetBx+2)*tan(2*pi/3);
H2C(:,offsetC+1) = H2C(:,offsetC+1) + H2C(:,offsetCx+3)*tan(4*pi/3);
H2C(:,no_col) = [];
H3A(:,offsetB+1) = H3A(:,offsetB+1) + H3A(:,offsetBx+2)*tan(2*pi/3);
H3A(:,offsetC+1) = H3A(:,offsetC+1) + H3A(:,offsetCx+3)*tan(4*pi/3);
H3A(:,no_col) = [];
H3B(:,offsetB+1) = H3B(:,offsetB+1) + H3B(:,offsetBx+2)*tan(2*pi/3);
H3B(:,offsetC+1) = H3B(:,offsetC+1) + H3B(:,offsetCx+3)*tan(4*pi/3);
H3B(:,no_col) = [];
H3C(:,offsetB+1) = H3C(:,offsetB+1) + H3C(:,offsetBx+2)*tan(2*pi/3);
H3C(:,offsetC+1) = H3C(:,offsetC+1) + H3C(:,offsetCx+3)*tan(4*pi/3);
H3C(:,no_col) = [];

% 4-5) Power Flow Measurements
H4A = zeros(npfA,num_state_var);
H5A = zeros(nqfA,num_state_var);
H4B = zeros(npfB,num_state_var);
H5B = zeros(nqfB,num_state_var);
H4C = zeros(npfC,num_state_var);
H5C = zeros(nqfC,num_state_var);
for i=1:npfA
    mA = fbuspfA(i);
    nA = tbuspfA(i);
    H4A(i,mA) = - G_state(mA,nA);
    H4A(i,nA) = G_state(mA,nA);
    H4A(i,mA + offsetAx) = B_state(mA,nA);
    H4A(i,nA + offsetAx) = - B_state(mA,nA);
    H5A(i,mA) = - B_state(mA,nA);
    H5A(i,nA) = B_state(mA,nA);
    H5A(i,mA + offsetAx) = - G_state(mA,nA);
    H5A(i,nA + offsetAx) = G_state(mA,nA);
    index = branch_cod_A==brpfA(i);
    if branch_cod_B(index) ~= 0
        mB = branch_first_B(index) + offsetB;
        nB = branch_end_B(index) + offsetB;
        H4A(i,mB) = - G_state(mA,nB);
        H4A(i,nB) = G_state(mA,nB);
        H4A(i,mB + offsetAx) = B_state(mA,nB);
        H4A(i,nB + offsetAx) = - B_state(mA,nB);
        H5A(i,mB) = - B_state(mA,nB);
        H5A(i,nB) = B_state(mA,nB);
        H5A(i,mB + offsetAx) = - G_state(mA,nB);
        H5A(i,nB + offsetAx) = G_state(mA,nB);
    end
    if branch_cod_C(index) ~= 0
        mC = branch_first_C(index) + offsetC;
        nC = branch_end_C(index) + offsetC;
        H4A(i,mC) = - G_state(mA,nC);
        H4A(i,nC) = G_state(mA,nC);
        H4A(i,mC + offsetAx) = B_state(mA,nC);
        H4A(i,nC + offsetAx) = - B_state(mA,nC);
        H5A(i,mC) = - B_state(mA,nC);
        H5A(i,nC) = B_state(mA,nC);
        H5A(i,mC + offsetAx) = - G_state(mA,nC);
        H5A(i,nC + offsetAx) = G_state(mA,nC);
    end
end
for i=1:npfB
    idx_p = pfstartB-1 + i;
    idx_q = qfstartB-1 + i;
    rot_mat = [cos(2*pi/3) , sin(2*pi/3); sin(2*pi/3), -cos(2*pi/3)];
    covar_pq_eq = rot_mat * [sigma2B(idx_p), 0; 0, sigma2B(idx_q)] * rot_mat';
    covariancesB(idx_p) = covar_pq_eq(1,1);
    covariancesB(idx_q) = covar_pq_eq(2,2);
    hidxB = [hidxB, idx_p, idx_q];
    vidxB = [vidxB, idx_q, idx_p];
    covariancesB = [covariancesB; covar_pq_eq(1,2); covar_pq_eq(2,1)];
    mB = fbuspfB(i) + offsetB;
    nB = tbuspfB(i) + offsetB;
    H4B(i,mB) = - G_state(mB,nB);
    H4B(i,nB) = G_state(mB,nB);
    H4B(i,mB + offsetAx) = B_state(mB,nB);
    H4B(i,nB + offsetAx) = - B_state(mB,nB);
    H5B(i,mB) = - B_state(mB,nB);
    H5B(i,nB) = B_state(mB,nB);
    H5B(i,mB + offsetAx) = - G_state(mB,nB);
    H5B(i,nB + offsetAx) = G_state(mB,nB);
    index = branch_cod_B==brpfB(i);
    if branch_cod_A(index) ~= 0
        mA = branch_first_A(index);
        nA = branch_end_A(index);
        H4B(i,mA) = - G_state(mB,nA);
        H4B(i,nA) = G_state(mB,nA);
        H4B(i,mA + offsetAx) = B_state(mB,nA);
        H4B(i,nA + offsetAx) = - B_state(mB,nA);
        H5B(i,mA) = - B_state(mB,nA);
        H5B(i,nA) = B_state(mB,nA);
        H5B(i,mA + offsetAx) = - G_state(mB,nA);
        H5B(i,nA + offsetAx) = G_state(mB,nA);
    end
    if branch_cod_C(index) ~= 0
        mC = branch_first_C(index) + offsetC;
        nC = branch_end_C(index) + offsetC;
        H4B(i,mC) = - G_state(mB,nC);
        H4B(i,nC) = G_state(mB,nC);
        H4B(i,mC + offsetAx) = B_state(mB,nC);
        H4B(i,nC + offsetAx) = - B_state(mB,nC);
        H5B(i,mC) = - B_state(mB,nC);
        H5B(i,nC) = B_state(mB,nC);
        H5B(i,mC + offsetAx) = - G_state(mB,nC);
        H5B(i,nC + offsetAx) = G_state(mB,nC);
    end
end
for i=1:npfC
    idx_p = pfstartC-1 + i;
    idx_q = qfstartC-1 + i;
    rot_mat = [cos(4*pi/3) , sin(4*pi/3); sin(4*pi/3), -cos(4*pi/3)];
    covar_pq_eq = rot_mat * [sigma2C(idx_p), 0; 0, sigma2C(idx_q)] * rot_mat';
    covariancesC(idx_p) = covar_pq_eq(1,1);
    covariancesC(idx_q) = covar_pq_eq(2,2);
    hidxC = [hidxC, idx_p, idx_q];
    vidxC = [vidxC, idx_q, idx_p];
    covariancesC = [covariancesC; covar_pq_eq(1,2); covar_pq_eq(2,1)];
    mC = fbuspfC(i) + offsetC;
    nC = tbuspfC(i) + offsetC;
    H4C(i,mC) = - G_state(mC,nC);
    H4C(i,nC) = G_state(mC,nC);
    H4C(i,mC + offsetAx) = B_state(mC,nC);
    H4C(i,nC + offsetAx) = - B_state(mC,nC);
    H5C(i,mC) = - B_state(mC,nC);
    H5C(i,nC) = B_state(mC,nC);
    H5C(i,mC + offsetAx) = - G_state(mC,nC);
    H5C(i,nC + offsetAx) = G_state(mC,nC);
    index = branch_cod_C==brpfC(i);
    if branch_cod_A(index) ~= 0
        mA = branch_first_A(index);
        nA = branch_end_A(index);
        H4C(i,mA) = - G_state(mC,nA);
        H4C(i,nA) = G_state(mC,nA);
        H4C(i,mA + offsetAx) = B_state(mC,nA);
        H4C(i,nA + offsetAx) = - B_state(mC,nA);
        H5C(i,mA) = - B_state(mC,nA);
        H5C(i,nA) = B_state(mC,nA);
        H5C(i,mA + offsetAx) = - G_state(mC,nA);
        H5C(i,nA + offsetAx) = G_state(mC,nA);
    end
    if branch_cod_B(index) ~= 0
        mB = branch_first_B(index) + offsetB;
        nB = branch_end_B(index) + offsetB;
        H4C(i,mB) = - G_state(mC,nB);
        H4C(i,nB) = G_state(mC,nB);
        H4C(i,mB + offsetAx) = B_state(mC,nB);
        H4C(i,nB + offsetAx) = - B_state(mC,nB);
        H5C(i,mB) = - B_state(mC,nB);
        H5C(i,nB) = B_state(mC,nB);
        H5C(i,mB + offsetAx) = - G_state(mC,nB);
        H5C(i,nB + offsetAx) = G_state(mC,nB);
    end
end
H4A(:,offsetB+1) = H4A(:,offsetB+1) + H4A(:,offsetBx+2)*tan(2*pi/3);
H4A(:,offsetC+1) = H4A(:,offsetC+1) + H4A(:,offsetCx+3)*tan(4*pi/3);
H5A(:,offsetB+1) = H5A(:,offsetB+1) + H5A(:,offsetBx+2)*tan(2*pi/3);
H5A(:,offsetC+1) = H5A(:,offsetC+1) + H5A(:,offsetCx+3)*tan(4*pi/3);
H4A(:,no_col) = [];
H5A(:,no_col) = [];
H4B(:,offsetB+1) = H4B(:,offsetB+1) + H4B(:,offsetBx+2)*tan(2*pi/3);
H4B(:,offsetC+1) = H4B(:,offsetC+1) + H4B(:,offsetCx+3)*tan(4*pi/3);
H5B(:,offsetB+1) = H5B(:,offsetB+1) + H5B(:,offsetBx+2)*tan(2*pi/3);
H5B(:,offsetC+1) = H5B(:,offsetC+1) + H5B(:,offsetCx+3)*tan(4*pi/3);
H4B(:,no_col) = [];
H5B(:,no_col) = [];
H4C(:,offsetB+1) = H4C(:,offsetB+1) + H4C(:,offsetBx+2)*tan(2*pi/3);
H4C(:,offsetC+1) = H4C(:,offsetC+1) + H4C(:,offsetCx+3)*tan(4*pi/3);
H5C(:,offsetB+1) = H5C(:,offsetB+1) + H5C(:,offsetBx+2)*tan(2*pi/3);
H5C(:,offsetC+1) = H5C(:,offsetC+1) + H5C(:,offsetCx+3)*tan(4*pi/3);
H4C(:,no_col) = [];
H5C(:,no_col) = [];

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

% Initial state Phase A,B,C
State = [VrA;VrB;VrC;VxA(2:end);VxB(2:end);VxC(2:end)];

num_iter = 0;
epsilonV = 5;
while (epsilonV > 10^-7 && num_iter <= 100)
    
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
    h1A = V_absA(busviA); 
    H1A = zeros(nviA,num_state_var);
    for i=1:nviA
        index = busviA(i);
        H1A(i, index) = cos(V_thetaA(index));
        H1A(i, index + offsetAx) = sin(V_thetaA(index));
    end
    H1A(:,no_col) = [];
    %Phase B
    h1B = V_absB(busviB); 
    H1B = zeros(nviB,num_state_var);
    for i=1:nviB
        index = busviB(i);
        H1B(i, index + offsetB) = cos(V_thetaB(index));
        H1B(i, index + offsetBx+1) = sin(V_thetaB(index));
    end
    H1B(:,offsetB+1) = H1B(:,offsetB+1) + H1B(:,offsetBx+2)*tan(2*pi/3);
    H1B(:,offsetC+1) = H1B(:,offsetC+1) + H1B(:,offsetCx+3)*tan(4*pi/3);
    H1B(:,no_col) = [];
    %Phase C
    h1C = V_absC(busviC);
    H1C = zeros(nviC,num_state_var);
    for i=1:nviC
        index = busviC(i);
        H1C(i, index + offsetC) = cos(V_thetaC(index));
        H1C(i, index + offsetCx+2) = sin(V_thetaC(index));
    end
    H1C(:,offsetB+1) = H1C(:,offsetB+1) + H1C(:,offsetBx+2)*tan(2*pi/3);
    H1C(:,offsetC+1) = H1C(:,offsetC+1) + H1C(:,offsetCx+3)*tan(4*pi/3);
    H1C(:,no_col) = [];
  
    %%% 2-3) Power Injection Measurements
    h2A = H2A*State;
    h2B = H2B*State;
    h2C = H2C*State; 
    h3A = H3A*State;
    h3B = H3B*State;
    h3C = H3C*State;
   
    %%% 4-5) Power Flow Measurements
    h4A = H4A*State;
    h4B = H4B*State;
    h4C = H4C*State;
    h5A = H5A*State;
    h5B = H5B*State;
    h5C = H5C*State;
    
    %%% 6) Current Magnitude Measurements
    % Phase A
    h6A_re = zeros(niampfA,1);
    h6A_im = zeros(niampfA,1);
    h6A = zeros(niampfA,1);
    H6A = zeros(niampfA,num_state_var);
    for i=1:niampfA
        index = branch_cod_A==briampfA(i);
        mA = fbusiampfA(i);
        nA = tbusiampfA(i);
        h6A_re(i,1) = Yabs_state(mA,nA)*((Vr_red(nA)-Vr_red(mA))*cos(Yphase_state(mA,nA)) + (Vx_red(mA)-Vx_red(nA))*sin(Yphase_state(mA,nA)));
        h6A_im(i,1) = Yabs_state(mA,nA)*((Vr_red(nA)-Vr_red(mA))*sin(Yphase_state(mA,nA)) + (Vx_red(nA)-Vx_red(mA))*cos(Yphase_state(mA,nA)));
        if branch_cod_B(index) ~= 0
            mB = branch_first_B(index) + offsetB;
            nB = branch_end_B(index) + offsetB;
            h6A_re(i,1) = h6A_re(i,1) + Yabs_state(mA,nB)*((Vr_red(nB)-Vr_red(mB))*cos(Yphase_state(mA,nB)) + (Vx_red(mB)-Vx_red(nB))*sin(Yphase_state(mA,nB)));
            h6A_im(i,1) = h6A_im(i,1) + Yabs_state(mA,nB)*((Vr_red(nB)-Vr_red(mB))*sin(Yphase_state(mA,nB)) + (Vx_red(nB)-Vx_red(mB))*cos(Yphase_state(mA,nB)));
        end
        if branch_cod_C(index) ~= 0
            mC = branch_first_C(index) + offsetC;
            nC = branch_end_C(index) + offsetC;
            h6A_re(i,1) = h6A_re(i,1) + Yabs_state(mA,nC)*((Vr_red(nC)-Vr_red(mC))*cos(Yphase_state(mA,nC)) + (Vx_red(mC)-Vx_red(nC))*sin(Yphase_state(mA,nC)));
            h6A_im(i,1) = h6A_im(i,1) + Yabs_state(mA,nC)*((Vr_red(nC)-Vr_red(mC))*sin(Yphase_state(mA,nC)) + (Vx_red(nC)-Vx_red(mC))*cos(Yphase_state(mA,nC)));
        end
        h6A(i,1) = sqrt(h6A_re(i,1)^2 + h6A_im(i,1)^2);
        H6A(i,mA) = - Yabs_state(mA,nA)*(cos(Yphase_state(mA,nA))*h6A_re(i,1)+sin(Yphase_state(mA,nA))*h6A_im(i,1))/h6A(i,1);
        H6A(i,nA) = - H6A(i,mA);
        H6A(i,mA + offsetAx) = - Yabs_state(mA,nA)*(cos(Yphase_state(mA,nA))*h6A_im(i,1)-sin(Yphase_state(mA,nA))*h6A_re(i,1))/h6A(i,1);
        H6A(i,nA + offsetAx) = - H6A(i,mA + offsetAx);
        if branch_cod_B(index) ~= 0
            H6A(i,mB) = - Yabs_state(mA,nB)*(cos(Yphase_state(mA,nB))*h6A_re(i,1)+sin(Yphase_state(mA,nB))*h6A_im(i,1))/h6A(i,1);
            H6A(i,nB) = - H6A(i,mB);
            H6A(i,mB + offsetAx) = - Yabs_state(mA,nB)*(cos(Yphase_state(mA,nB))*h6A_im(i,1)-sin(Yphase_state(mA,nB))*h6A_re(i,1))/h6A(i,1);
            H6A(i,nB + offsetAx) = - H6A(i,mB + offsetAx);
        end
        if branch_cod_C(index) ~= 0
            H6A(i,mC) = - Yabs_state(mA,nC)*(cos(Yphase_state(mA,nC))*h6A_re(i,1)+sin(Yphase_state(mA,nC))*h6A_im(i,1))/h6A(i,1);
            H6A(i,nC) = - H6A(i,mC);
            H6A(i,mC + offsetAx) = - Yabs_state(mA,nC)*(cos(Yphase_state(mA,nC))*h6A_im(i,1)-sin(Yphase_state(mA,nC))*h6A_re(i,1))/h6A(i,1);
            H6A(i,nC + offsetAx) = - H6A(i,mC + offsetAx);
        end
    end
    H6A(:,offsetB+1) = H6A(:,offsetB+1) + H6A(:,offsetBx+2)*tan(2*pi/3);
    H6A(:,offsetC+1) = H6A(:,offsetC+1) + H6A(:,offsetCx+3)*tan(4*pi/3);
    H6A(:,no_col) = [];
    % Phase B
    h6B_re = zeros(niampfB,1);
    h6B_im = zeros(niampfB,1);
    h6B = zeros(niampfB,1);
    H6B = zeros(niampfB,num_state_var);
    for i=1:niampfB
        index = branch_cod_B==briampfB(i);
        mB = fbusiampfB(i) + offsetB;
        nB = tbusiampfB(i) + offsetB;
        h6B_re(i,1) = Yabs_state(mB,nB)*((Vr_red(nB)-Vr_red(mB))*cos(Yphase_state(mB,nB)) + (Vx_red(mB)-Vx_red(nB))*sin(Yphase_state(mB,nB)));
        h6B_im(i,1) = Yabs_state(mB,nB)*((Vr_red(nB)-Vr_red(mB))*sin(Yphase_state(mB,nB)) + (Vx_red(nB)-Vx_red(mB))*cos(Yphase_state(mB,nB)));
        if branch_cod_A(index) ~= 0
            mA = branch_first_A(index);
            nA = branch_end_A(index);
            h6B_re(i,1) = h6B_re(i,1) + Yabs_state(mB,nA)*((Vr_red(nA)-Vr_red(mA))*cos(Yphase_state(mB,nA)) + (Vx_red(mA)-Vx_red(nA))*sin(Yphase_state(mB,nA)));
            h6B_im(i,1) = h6B_im(i,1) + Yabs_state(mB,nA)*((Vr_red(nA)-Vr_red(mA))*sin(Yphase_state(mB,nA)) + (Vx_red(nA)-Vx_red(mA))*cos(Yphase_state(mB,nA)));
        end
        if branch_cod_C(index) ~= 0
            mC = branch_first_C(index) + offsetC;
            nC = branch_end_C(index) + offsetC;
            h6B_re(i,1) = h6B_re(i,1) + Yabs_state(mB,nC)*((Vr_red(nC)-Vr_red(mC))*cos(Yphase_state(mB,nC)) + (Vx_red(mC)-Vx_red(nC))*sin(Yphase_state(mB,nC)));
            h6B_im(i,1) = h6B_im(i,1) + Yabs_state(mB,nC)*((Vr_red(nC)-Vr_red(mC))*sin(Yphase_state(mB,nC)) + (Vx_red(nC)-Vx_red(mC))*cos(Yphase_state(mB,nC)));
        end
        h6B(i,1) = sqrt(h6B_re(i,1)^2 + h6B_im(i,1)^2);
        H6B(i,mB) = - Yabs_state(mB,nB)*(cos(Yphase_state(mB,nB))*h6B_re(i,1)+sin(Yphase_state(mB,nB))*h6B_im(i,1))/h6B(i,1);
        H6B(i,nB) = - H6B(i,mB);
        H6B(i,mB + offsetAx) = - Yabs_state(mB,nB)*(cos(Yphase_state(mB,nB))*h6B_im(i,1)-sin(Yphase_state(mB,nB))*h6B_re(i,1))/h6B(i,1);
        H6B(i,nB + offsetAx) = - H6B(i,mB + offsetAx);
        if branch_cod_A(index) ~= 0
            H6B(i,mA) = - Yabs_state(mB,nA)*(cos(Yphase_state(mB,nA))*h6B_re(i,1)+sin(Yphase_state(mB,nA))*h6B_im(i,1))/h6B(i,1);
            H6B(i,nA) = - H6B(i,mA);
            H6B(i,mA + offsetAx) = - Yabs_state(mB,nA)*(cos(Yphase_state(mB,nA))*h6B_im(i,1)-sin(Yphase_state(mB,nA))*h6B_re(i,1))/h6B(i,1);
            H6B(i,nA + offsetAx) = - H6B(i,mA + offsetAx);
        end
        if branch_cod_C(index) ~= 0
            H6B(i,mC) = - Yabs_state(mB,nC)*(cos(Yphase_state(mB,nC))*h6B_re(i,1)+sin(Yphase_state(mB,nC))*h6B_im(i,1))/h6B(i,1);
            H6B(i,nC) = - H6B(i,mC);
            H6B(i,mC + offsetAx) = - Yabs_state(mB,nC)*(cos(Yphase_state(mB,nC))*h6B_im(i,1)-sin(Yphase_state(mB,nC))*h6B_re(i,1))/h6B(i,1);
            H6B(i,nC + offsetAx) = - H6B(i,mC + offsetAx);
        end
    end
    H6B(:,offsetB+1) = H6B(:,offsetB+1) + H6B(:,offsetBx+2)*tan(2*pi/3);
    H6B(:,offsetC+1) = H6B(:,offsetC+1) + H6B(:,offsetCx+3)*tan(4*pi/3);
    H6B(:,no_col) = [];
    % Phase C
    h6C_re = zeros(niampfC,1);
    h6C_im = zeros(niampfC,1);
    h6C = zeros(niampfC,1);
    H6C = zeros(niampfC,num_state_var);
    for i=1:niampfC
        index = branch_cod_C==briampfC(i);
        mC = fbusiampfC(i) + offsetC;
        nC = tbusiampfC(i) + offsetC;
        h6C_re(i,1) = Yabs_state(mC,nC)*((Vr_red(nC)-Vr_red(mC))*cos(Yphase_state(mC,nC)) + (Vx_red(mC)-Vx_red(nC))*sin(Yphase_state(mC,nC)));
        h6C_im(i,1) = Yabs_state(mC,nC)*((Vr_red(nC)-Vr_red(mC))*sin(Yphase_state(mC,nC)) + (Vx_red(nC)-Vx_red(mC))*cos(Yphase_state(mC,nC)));
        if branch_cod_A(index) ~= 0
            mA = branch_first_A(index);
            nA = branch_end_A(index);
            h6C_re(i,1) = h6C_re(i,1) + Yabs_state(mC,nA)*((Vr_red(nA)-Vr_red(mA))*cos(Yphase_state(mC,nA)) + (Vx_red(mA)-Vx_red(nA))*sin(Yphase_state(mC,nA)));
            h6C_im(i,1) = h6C_im(i,1) + Yabs_state(mC,nA)*((Vr_red(nA)-Vr_red(mA))*sin(Yphase_state(mC,nA)) + (Vx_red(nA)-Vx_red(mA))*cos(Yphase_state(mC,nA)));
        end
        if branch_cod_B(index) ~= 0
            mB = branch_first_B(index) + offsetB;
            nB = branch_end_B(index) + offsetB;
            h6C_re(i,1) = h6C_re(i,1) + Yabs_state(mC,nB)*((Vr_red(nB)-Vr_red(mB))*cos(Yphase_state(mC,nB)) + (Vx_red(mB)-Vx_red(nB))*sin(Yphase_state(mC,nB)));
            h6C_im(i,1) = h6C_im(i,1) + Yabs_state(mC,nB)*((Vr_red(nB)-Vr_red(mB))*sin(Yphase_state(mC,nB)) + (Vx_red(nB)-Vx_red(mB))*cos(Yphase_state(mC,nB)));
        end
        h6C(i,1) = sqrt(h6C_re(i,1)^2 + h6C_im(i,1)^2);
        H6C(i,mC) = - Yabs_state(mC,nC)*(cos(Yphase_state(mC,nC))*h6C_re(i,1)+sin(Yphase_state(mC,nC))*h6C_im(i,1))/h6C(i,1);
        H6C(i,nC) = - H6C(i,mC);
        H6C(i,mC + offsetAx) = - Yabs_state(mC,nC)*(cos(Yphase_state(mC,nC))*h6C_im(i,1)-sin(Yphase_state(mC,nC))*h6C_re(i,1))/h6C(i,1);
        H6C(i,nC + offsetAx) = - H6C(i,mC + offsetAx);
        if branch_cod_A(index) ~= 0
            H6C(i,mA) = - Yabs_state(mC,nA)*(cos(Yphase_state(mC,nA))*h6C_re(i,1)+sin(Yphase_state(mC,nA))*h6C_im(i,1))/h6C(i,1);
            H6C(i,nA) = - H6C(i,mA);
            H6C(i,mA + offsetAx) = - Yabs_state(mC,nA)*(cos(Yphase_state(mC,nA))*h6C_im(i,1)-sin(Yphase_state(mC,nA))*h6C_re(i,1))/h6C(i,1);
            H6C(i,nA + offsetAx) = - H6C(i,mA + offsetAx);
        end
        if branch_cod_B(index) ~= 0
            H6C(i,mB) = - Yabs_state(mC,nB)*(cos(Yphase_state(mC,nB))*h6C_re(i,1)+sin(Yphase_state(mC,nB))*h6C_im(i,1))/h6C(i,1);
            H6C(i,nB) = - H6C(i,mB);
            H6C(i,mB + offsetAx) = - Yabs_state(mC,nB)*(cos(Yphase_state(mC,nB))*h6C_im(i,1)-sin(Yphase_state(mC,nB))*h6C_re(i,1))/h6C(i,1);
            H6C(i,nB + offsetAx) = - H6C(i,mB + offsetAx);
        end
    end
    H6C(:,offsetB+1) = H6C(:,offsetB+1) + H6C(:,offsetBx+2)*tan(2*pi/3);
    H6C(:,offsetC+1) = H6C(:,offsetC+1) + H6C(:,offsetCx+3)*tan(4*pi/3);
    H6C(:,no_col) = []; 
        
    %%%  fine misure
    
    % Composizione vettore misure e stima
    % Phase A
    hA = [h1A; h2A; h3A; h4A; h5A; h6A];
    HA = [H1A; H2A; H3A; H4A; H5A; H6A];
     % Phase B
    hB = [h1B; h2B; h3B; h4B; h5B; h6B];
    HB = [H1B; H2B; H3B; H4B; H5B; H6B];
     % Phase C
    hC = [h1C; h2C; h3C; h4C; h5C; h6C];
    HC = [H1C; H2C; H3C; H4C; H5C; H6C];
    
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
    VrA = State(1:offsetB);
    VxA(2:end) = State(offsetAx+1:offsetBx);
    % Phase B
    VrB = State(offsetB+1:offsetC);
    VxB(2:end) = State(offsetBx+1:offsetCx);
    VxB(1) = VrB(1)*tan(2*pi/3);
    % Phase C
    VrC = State(offsetC+1:offsetAx);
    VxC(2:end) = State(offsetCx+1:end);
    VxC(1) = VrC(1)*tan(4*pi/3);
    
    Vr_red = [VrA; VrB; VrC];
    Vx_red = [VxA; VxB; VxC];
    
    VA = complex(VrA,VxA);
    VB = complex(VrB,VxB);
    VC = complex(VrC,VxC);
    
    V_absA = abs(VA);
    V_absB = abs(VB);
    V_absC = abs(VC);
    
    V_thetaA = angle(VA);
    V_thetaB = angle(VB);
    V_thetaC = angle(VC);    

    num_iter = num_iter + 1;
    epsilonV = max(abs(Delta_X));

end

% Vincenzo
Vr = zeros(size(Vr));
Vx = zeros(size(Vx));

Vr(A_nodes_indices) = VrA;
Vr(B_nodes_indices+num_nodes) = VrB;
Vr(C_nodes_indices+2*num_nodes) = VrC;

Vx(A_nodes_indices) = VxA;
Vx(B_nodes_indices+num_nodes) = VxB;
Vx(C_nodes_indices+2*num_nodes) = VxC;

V = complex(Vr,Vx);

I = zeros(3*num_branches,1);
for i = 1:num_branches
    i2 = i+num_branches;
    i3 = i+2*num_branches;
    m = branch_first(i);
    n = branch_end(i);
    I(i,1) = 0;
    I(i2,1) = 0;
    I(i3,1) = 0;
    for j = 1:3
        k = (j-1)*num_nodes;
        I(i,1) = I(i,1) -(V(m+k) - V(n+k))*(Y_tot(m,n+k));
        I(i2,1) = I(i2,1) -(V(m+k) - V(n+k))*(Y_tot(m+num_nodes,n+k));
        I(i3,1) = I(i3,1) -(V(m+k) - V(n+k))*(Y_tot(m+2*num_nodes,n+k));
    end
end
Ir = real(I);
Ix = imag(I);
