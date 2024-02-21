function[V_true,I_br_true,S_br_true,S_inj_true,num_iter]=Power_flow_New_Raph_3Ph(num_nodes, node_3ph, node_3ph_reverse, Y, pseudomeas, nominal_linedata, Zbr_set, V, theta)

%%%  Il vettore P deve contenere le P note relative ai nodi PV e PQ mentre
%%%  il vettore Q deve contenere le Q note relative ai soli nodi PQ.
%%%  Le potenze assorbite ai nodi (carico) vanno indicate con segno
%%%  negativo, mentre quelle iniettate (generazione) con segno positivo.
%%%  Sia i vettori P che Q che V che theta devono essere vettori colonna.
%%%  Nella definizione dei nodi (e quindi dei vettori V e theta), mettere
%%%  all'indice 1 il nodo di slack, seguito da tutti i nodi PQ e infine dai nodi PV.

branch_first = nominal_linedata(:, 1);
branch_first_A = nominal_linedata(:,3);
branch_first_B = nominal_linedata(:,4);
branch_first_C = nominal_linedata(:,5);
branch_end_A = nominal_linedata(:,6);
branch_end_B = nominal_linedata(:,7);
branch_end_C = nominal_linedata(:,8);
% branch_cod_A = nominal_linedata(:,10);
% branch_cod_B = nominal_linedata(:,11);
% branch_cod_C = nominal_linedata(:,12);
% num_br_A = max(branch_cod_A);
% num_br_B = max(branch_cod_B);
% num_br_C = max(branch_cod_C);

num_A_nodes = max(branch_end_A);
num_B_nodes = max(branch_end_B);
num_C_nodes = max(branch_end_C);

existing_nodes  = (node_3ph ~= 0);
indices = 1:num_nodes;
VA = V(indices);
VB = V(indices + num_nodes);
VC = V(indices + 2*num_nodes);
thetaA = theta(indices);
thetaB = theta(indices + num_nodes);
thetaC = theta(indices + 2*num_nodes);
V = [VA(existing_nodes(:,1));VB(existing_nodes(:,2));VC(existing_nodes(:,3))];
theta = [thetaA(existing_nodes(:,1));thetaB(existing_nodes(:,2));thetaC(existing_nodes(:,3))];

num_state_nodes = num_A_nodes + num_B_nodes + num_C_nodes;
% num_branches = num_br_A + num_br_B + num_br_C;
num_branches = length(branch_first);

% nodes_with_A = node_3ph(node_3ph(:,1) > 0, 1); % indici dei nodi che hanno fase A nel vettore della fase A
% nodes_with_B = node_3ph(node_3ph(:,2) > 0, 2);
% nodes_with_C = node_3ph(node_3ph(:,3) > 0, 3);


% nodi della rete che hanno le rispettive fasi
nodes_with_A = node_3ph_reverse(1:num_A_nodes,1);
nodes_with_A = nodes_with_A(2:end); % escludo il primo che non riguardera le potenze
nodes_with_B = node_3ph_reverse(1:num_B_nodes,2);
nodes_with_B = nodes_with_B(2:end); % escludo il primo che non riguardera le potenze
nodes_with_C = node_3ph_reverse(1:num_C_nodes,3);
nodes_with_C = nodes_with_C(2:end); % escludo il primo che non riguardera le potenze

PA = pseudomeas(nodes_with_A-1,1); % tolgo 1 perche le pseudo si riferiscono ai nodi dal 2 in poi
PB = pseudomeas(nodes_with_B-1,2);
PC = pseudomeas(nodes_with_C-1,3);
QA = pseudomeas(nodes_with_A-1,4);
QB = pseudomeas(nodes_with_B-1,5);
QC = pseudomeas(nodes_with_C-1,6);

G=real(Y);
B=imag(Y);
S=[PA;QA;PB;QB;PC;QC];
epsilon=1;
num_iter=0;
Delta_V = zeros(num_state_nodes,1);
while epsilon > 10^-6               %%% criterio convergenza %%%
    D_P=zeros(num_nodes-1,3);
    D_Q=zeros(num_nodes-1,3);
    
    for i=2:num_nodes
        A_n_i = node_3ph(i,1);
        B_n_i = node_3ph(i,2) + num_A_nodes;
        C_n_i = node_3ph(i,3) + num_A_nodes + num_B_nodes;
        if (A_n_i > 0)
            D_P(i-1,1)=V(A_n_i)*(G(A_n_i,:)*(V(:).*cos(theta(A_n_i)-theta))+B(A_n_i,:)*(V(:).*sin(theta(A_n_i)-theta)));
            D_Q(i-1,1)=V(A_n_i)*(G(A_n_i,:)*(V(:).*sin(theta(A_n_i)-theta))-B(A_n_i,:)*(V(:).*cos(theta(A_n_i)-theta)));
        end
        if (B_n_i > 0)
            D_P(i-1,2)=V(B_n_i)*(G(B_n_i,:)*(V(:).*cos(theta(B_n_i)-theta))+B(B_n_i,:)*(V(:).*sin(theta(B_n_i)-theta)));
            D_Q(i-1,2)=V(B_n_i)*(G(B_n_i,:)*(V(:).*sin(theta(B_n_i)-theta))-B(B_n_i,:)*(V(:).*cos(theta(B_n_i)-theta)));
        end
        if (C_n_i > 0)
            D_P(i-1,3)=V(C_n_i)*(G(C_n_i,:)*(V(:).*cos(theta(C_n_i)-theta))+B(C_n_i,:)*(V(:).*sin(theta(C_n_i)-theta)));
            D_Q(i-1,3)=V(C_n_i)*(G(C_n_i,:)*(V(:).*sin(theta(C_n_i)-theta))-B(C_n_i,:)*(V(:).*cos(theta(C_n_i)-theta)));
        end

%         for j=1:num_state_nodes  
%             if (A_n_i > 0)
%                 D_P(i-1,1)=D_P(i-1,1) + V(A_n_i)*V(j)*(G(A_n_i,j)*cos(theta(A_n_i)-theta(j))+B(A_n_i,j)*sin(theta(A_n_i)-theta(j)));
%                 D_Q(i-1,1)=D_Q(i-1,1) + V(A_n_i)*V(j)*(G(A_n_i,j)*sin(theta(A_n_i)-theta(j))-B(A_n_i,j)*cos(theta(A_n_i)-theta(j)));
%             end
%             if (B_n_i - num_A_nodes > 0)
%                 D_P(i-1,2)=D_P(i-1,2) + V(B_n_i)*V(j)*(G(B_n_i,j)*cos(theta(B_n_i)-theta(j))+B(B_n_i,j)*sin(theta(B_n_i)-theta(j)));
%                 D_Q(i-1,2)=D_Q(i-1,2) + V(B_n_i)*V(j)*(G(B_n_i,j)*sin(theta(B_n_i)-theta(j))-B(B_n_i,j)*cos(theta(B_n_i)-theta(j)));
%             end
%             if (C_n_i - num_A_nodes - num_B_nodes > 0)
%                 D_P(i-1,3)=D_P(i-1,3) + V(C_n_i)*V(j)*(G(C_n_i,j)*cos(theta(C_n_i)-theta(j))+B(C_n_i,j)*sin(theta(C_n_i)-theta(j)));
%                 D_Q(i-1,3)=D_Q(i-1,3) + V(C_n_i)*V(j)*(G(C_n_i,j)*sin(theta(C_n_i)-theta(j))-B(C_n_i,j)*cos(theta(C_n_i)-theta(j)));
%             end
%         end
        %pippo=1;
    end
    
    D_S=[D_P(nodes_with_A-1,1); D_Q(nodes_with_A-1,1); D_P(nodes_with_B-1,2); D_Q(nodes_with_B-1,2); D_P(nodes_with_C-1,3); D_Q(nodes_with_C-1,3)];
    Delta_S = S-D_S;               %%% residuo della funzione obiettivo %%%
    
    epsilon = max(abs(Delta_S));
    
    num_iter = num_iter + 1;
    
         
    for i = nodes_with_A'          %%% derivate di P e Q rispetto a theta %%%
        A_n_i = node_3ph(i,1);    % indice del nodo nel settore dello stato relativo alla fase A 

        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A
               state_idx = state_idx - 1;
               continue;
            end
            
            if A_n_i == j
                J11A(A_n_i - 1, state_idx) = -D_Q(i - 1, 1) - B(j,j)*V(j)^2;
                J21A(A_n_i - 1, state_idx) =  D_P(i - 1, 1) - G(j,j)*V(j)^2;
            else
                J11A(A_n_i - 1, state_idx)= V(A_n_i)*V(j)*(G(A_n_i,j)*sin(theta(A_n_i)-theta(j))-B(A_n_i,j)*cos(theta(A_n_i)-theta(j)));
                J21A(A_n_i - 1, state_idx)=-V(A_n_i)*V(j)*(G(A_n_i,j)*cos(theta(A_n_i)-theta(j))+B(A_n_i,j)*sin(theta(A_n_i)-theta(j)));
            end
        end
    end
    for i = nodes_with_B'          %%% derivate di PB e QB rispetto a theta %%%
        B_n_i = node_3ph(i,2) + num_A_nodes;    % indice del nodo nel settore dello stato relativo alla fase A 

        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
               state_idx = state_idx - 1;
               continue;
            end
            
            if B_n_i == j
                J11B(B_n_i - num_A_nodes - 1, state_idx) = -D_Q(i - 1, 2) - B(j,j)*V(j)^2;
                J21B(B_n_i - num_A_nodes - 1, state_idx) =  D_P(i - 1, 2) - G(j,j)*V(j)^2;
            else
                J11B(B_n_i - num_A_nodes - 1, state_idx)=  V(B_n_i)*V(j)*(G(B_n_i,j)*sin(theta(B_n_i)-theta(j))-B(B_n_i,j)*cos(theta(B_n_i)-theta(j)));
                J21B(B_n_i - num_A_nodes - 1, state_idx)= -V(B_n_i)*V(j)*(G(B_n_i,j)*cos(theta(B_n_i)-theta(j))+B(B_n_i,j)*sin(theta(B_n_i)-theta(j)));
            end
        end
    end
    for i = nodes_with_C'          %%% derivate di P e Q rispetto a theta %%%
        C_n_i = node_3ph(i,3) + num_A_nodes + num_B_nodes;    % indice del nodo nel settore dello stato relativo alla fase A 

        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
               state_idx = state_idx - 1;
               continue;
            end
            
            if C_n_i == j
                J11C(C_n_i - num_A_nodes - num_B_nodes - 1, state_idx) = -D_Q(i - 1, 3) - B(j,j)*V(j)^2;
                J21C(C_n_i - num_A_nodes - num_B_nodes - 1, state_idx) =  D_P(i - 1, 3) - G(j,j)*V(j)^2;
            else
                J11C(C_n_i - num_A_nodes - num_B_nodes - 1, state_idx)=  V(C_n_i)*V(j)*(G(C_n_i,j)*sin(theta(C_n_i)-theta(j))-B(C_n_i,j)*cos(theta(C_n_i)-theta(j)));
                J21C(C_n_i - num_A_nodes - num_B_nodes - 1, state_idx)= -V(C_n_i)*V(j)*(G(C_n_i,j)*cos(theta(C_n_i)-theta(j))+B(C_n_i,j)*sin(theta(C_n_i)-theta(j)));
            end
        end
    end
    
    
    for i=nodes_with_A'     %%% derivate di PA e QA rispetto a V %%%
        A_n_i = node_3ph(i,1);
        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
               state_idx = state_idx - 1;
               continue;
            end
            if A_n_i == j
                J12A(A_n_i-1,state_idx)=D_P(A_n_i-1,1) + G(j,j)*V(j)^2;
                J22A(A_n_i-1,state_idx)=D_Q(A_n_i-1,1) - B(j,j)*V(j)^2;
            else
                J12A(A_n_i-1,state_idx)=V(A_n_i)*V(j)*(G(A_n_i,j)*cos(theta(A_n_i)-theta(j))+B(A_n_i,j)*sin(theta(A_n_i)-theta(j)));
                J22A(A_n_i-1,state_idx)=V(A_n_i)*V(j)*(G(A_n_i,j)*sin(theta(A_n_i)-theta(j))-B(A_n_i,j)*cos(theta(A_n_i)-theta(j)));
            end
        end
    end
    for i=nodes_with_B'    %%% derivate di PB e QB rispetto a V %%%
        B_n_i = node_3ph(i,2) + num_A_nodes;
        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
               state_idx = state_idx - 1;
               continue;
            end
            if B_n_i == j
                J12B(B_n_i - num_A_nodes -1,state_idx) = D_P(B_n_i - num_A_nodes - 1,2) + G(j,j)*V(j)^2;
                J22B(B_n_i - num_A_nodes -1,state_idx) = D_Q(B_n_i - num_A_nodes - 1,2) - B(j,j)*V(j)^2;
            else
                J12B(B_n_i - num_A_nodes -1,state_idx) = V(B_n_i)*V(j)*(G(B_n_i,j)*cos(theta(B_n_i)-theta(j))+B(B_n_i,j)*sin(theta(B_n_i)-theta(j)));
                J22B(B_n_i - num_A_nodes -1,state_idx) = V(B_n_i)*V(j)*(G(B_n_i,j)*sin(theta(B_n_i)-theta(j))-B(B_n_i,j)*cos(theta(B_n_i)-theta(j)));
            end
        end
    end
    for i=nodes_with_C'     %%% derivate di PC e QC rispetto a V %%%
        C_n_i = node_3ph(i,3) + num_A_nodes + num_B_nodes;
        state_idx = 0;
        for j= 1 : num_state_nodes
            state_idx = state_idx + 1;
            if (j == 1 || j == num_A_nodes + 1 || j == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
               state_idx = state_idx - 1;
               continue;
            end
            if C_n_i == j
                J12C(C_n_i - num_A_nodes - num_B_nodes - 1,state_idx) = D_P(C_n_i - num_A_nodes - num_B_nodes - 1,3) + G(j,j)*V(j)^2;
                J22C(C_n_i - num_A_nodes - num_B_nodes - 1,state_idx) = D_Q(C_n_i - num_A_nodes - num_B_nodes - 1,3) - B(j,j)*V(j)^2;
            else
                J12C(C_n_i - num_A_nodes - num_B_nodes - 1,state_idx) = V(C_n_i)*V(j)*(G(C_n_i,j)*cos(theta(C_n_i)-theta(j))+B(C_n_i,j)*sin(theta(C_n_i)-theta(j)));
                J22C(C_n_i - num_A_nodes - num_B_nodes - 1,state_idx) = V(C_n_i)*V(j)*(G(C_n_i,j)*sin(theta(C_n_i)-theta(j))-B(C_n_i,j)*cos(theta(C_n_i)-theta(j)));
            end
        end
    end
    
    
    J=[J11A, J12A; J21A, J22A; J11B, J12B; J21B, J22B; J11C, J12C; J21C, J22C];
    Delta_X = J\Delta_S;
    
    state_idx = 0;
    
    for i = 1 : num_state_nodes
        state_idx = state_idx + 1;
        if (i == 1 || i == num_A_nodes + 1 || i == num_A_nodes + num_B_nodes + 1) % fase primo nodo A primo nodo B e primo nodo C
           state_idx = state_idx - 1;
           continue;
        end
        theta(i) = theta(i) + Delta_X(state_idx);
        Delta_V(i) = Delta_X(state_idx + num_state_nodes - 3);
    end
    V = V + Delta_V.*V;
end
    

V_abs_true = V;
V_theta_true = theta;
V_state_true = V_abs_true.*exp(1i*V_theta_true);

P_inj_true = zeros(num_nodes,3);
Q_inj_true = zeros(num_nodes,3);
for i=2:num_nodes
    A_n_i = node_3ph(i,1);
    B_n_i = node_3ph(i,2) + num_A_nodes;
    C_n_i = node_3ph(i,3) + num_A_nodes + num_B_nodes;
  
    if (A_n_i > 0)
        P_inj_true(i,1) = V(A_n_i)*(G(A_n_i,:)*(V(:).*cos(theta(A_n_i)-theta))+B(A_n_i,:)*(V(:).*sin(theta(A_n_i)-theta)));
        Q_inj_true(i,1) = V(A_n_i)*(G(A_n_i,:)*(V(:).*sin(theta(A_n_i)-theta))-B(A_n_i,:)*(V(:).*cos(theta(A_n_i)-theta)));
    end
    if (B_n_i - num_A_nodes > 0)
        P_inj_true(i,2) = V(B_n_i)*(G(B_n_i,:)*(V(:).*cos(theta(B_n_i)-theta))+B(B_n_i,:)*(V(:).*sin(theta(B_n_i)-theta)));
        Q_inj_true(i,2) = V(B_n_i)*(G(B_n_i,:)*(V(:).*sin(theta(B_n_i)-theta))-B(B_n_i,:)*(V(:).*cos(theta(B_n_i)-theta)));
    end
    if (C_n_i - num_A_nodes - num_B_nodes > 0)
        P_inj_true(i,3) = V(C_n_i)*(G(C_n_i,:)*(V(:).*cos(theta(C_n_i)-theta))+B(C_n_i,:)*(V(:).*sin(theta(C_n_i)-theta)));
        Q_inj_true(i,3) = V(C_n_i)*(G(C_n_i,:)*(V(:).*sin(theta(C_n_i)-theta))-B(C_n_i,:)*(V(:).*cos(theta(C_n_i)-theta)));
    end
end

S_inj_true = complex(P_inj_true,Q_inj_true);

isA = nominal_linedata(:,13);
isB = nominal_linedata(:,14);
isC = nominal_linedata(:,15);
I_br_true = zeros(num_branches, 3);
S_br_true = zeros(num_branches, 3);
for br_idx = 1:num_branches
    Z3true = Zbr_set{br_idx};
    iA = branch_first_A(br_idx);
    jA = branch_end_A(br_idx);
    iB = branch_first_B(br_idx) + num_A_nodes;
    jB = branch_end_B(br_idx) + num_A_nodes;
    iC = branch_first_C(br_idx) + num_A_nodes + num_B_nodes;
    jC = branch_end_C(br_idx) + num_A_nodes + num_B_nodes;
    
    Vfrom_br = [];
    Vto_br = [];
    columns = [];
    remove_index = 3;
    if (~isA(br_idx))
        Z3true(remove_index - 2, :) = [];
        Z3true(:, remove_index - 2) = [];
        remove_index = remove_index - 1;
    else
        Vfrom_br = [Vfrom_br; V_state_true(iA)];
        Vto_br = [Vto_br; V_state_true(jA)];
        columns = [columns, 1];
    end
    if (~isB(br_idx))
        Z3true(remove_index - 1, :) = [];
        Z3true(:, remove_index - 1) = [];
        remove_index = remove_index - 1;
    else
        Vfrom_br = [Vfrom_br; V_state_true(iB)];
        Vto_br = [Vto_br; V_state_true(jB)];
        columns = [columns, 2];
    end
    if (~isC(br_idx))
        Z3true(remove_index, :) = [];
        Z3true(:, remove_index) = [];
    else
        Vfrom_br = [Vfrom_br; V_state_true(iC)];
        Vto_br = [Vto_br; V_state_true(jC)];
        columns = [columns, 3];
    end
    
    I_br_true(br_idx, columns) = (Z3true\(Vfrom_br-Vto_br)).';
    S_br_true(br_idx, columns) = (Vfrom_br.').*conj(I_br_true(br_idx, columns));
end

V_true = zeros(num_nodes, 3);
for i=1:num_nodes
    A_n_i = node_3ph(i,1);
    B_n_i = node_3ph(i,2) + num_A_nodes;
    C_n_i = node_3ph(i,3) + num_A_nodes + num_B_nodes;
    if (A_n_i > 0)
        V_true(i,1) = V_state_true(A_n_i);
    end
    if (B_n_i - num_A_nodes > 0)
        V_true(i,2) = V_state_true(B_n_i);
    end
    if (C_n_i - num_A_nodes - num_B_nodes > 0)
        V_true(i,3) = V_state_true(C_n_i);
    end
end
