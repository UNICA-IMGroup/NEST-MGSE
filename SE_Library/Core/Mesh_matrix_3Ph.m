function [Gz, Gz_state, Gzmod, Gzangle, Gzmod_state, Gzangle_state, Gzcomplex, eq_curr_matrix] = Mesh_matrix_3Ph(num_nodes, nominal_linedata, num_branches, Zbr_set, zero_inj_cell, eq_curr_matrix)

branch_first = nominal_linedata(:,1);
branch_end = nominal_linedata(:,2);
branch_cod = nominal_linedata(:,9);
zeroinjA = zero_inj_cell{1,1}; zeroinjB = zero_inj_cell{1,2}; zeroinjC = zero_inj_cell{1,3};
nzeroinjA = length(zeroinjA); nzeroinjB = length(zeroinjB); nzeroinjC = length(zeroinjC);
branch_cod_A = nominal_linedata(:,10);
branch_cod_B = nominal_linedata(:,11);
branch_cod_C = nominal_linedata(:,12);
num_branchesA = sum(branch_cod_A > 0);
num_branchesB = sum(branch_cod_B > 0);
num_branchesC = sum(branch_cod_C > 0);
num_branchesA_red = num_branchesA - nzeroinjA;
num_branchesB_red = num_branchesB - nzeroinjB;
num_branchesC_red = num_branchesC - nzeroinjC;
eq_curr_matrix = full(eq_curr_matrix);


for i = 1:num_nodes
    a(:,i) = (branch_end == i);        %%% ogni colonna di a contiene un valore logico pari a 1 nella riga corrispondente al ramo che finisce sul nodo i
end
for i=1:num_nodes
    num_br_to_nod(i,1) = sum(a(:,i));   %%% il vettore risultante contiene il numero di rami che terminano sul nodo di indice i
end

nod_idx_mesh = find(num_br_to_nod>1);   %%% indice dei nodi che hanno piú di un ramo che converge su di essi
path_best_mesh = [];

for j=1:length(nod_idx_mesh)
    br_from = branch_first(branch_end == nod_idx_mesh(j));   %%% nodi di partenza dei rami che convergono sullo stesso nodo
    br_to_nod = br_from;
    var = 0;
    path_mesh = br_from;
    while var == 0                               
        var = 1;
        br_var2 = [];
        path_mesh2 =[];
        for i=1:length(br_from)
            if br_from(i) <= 1
                br_var = 0;
            else
                br_var = branch_first(branch_end == br_from(i));
                var = 0;
            end
            br_var2 = [br_var2; br_var];
            path_mesh2 = [path_mesh2; ones(length(br_var),1)*path_mesh(i,:)]; 
        end
        br_from = br_var2;
        path_mesh = [path_mesh2, br_from];    %%% formazione dei percorsi relativi alla maglia j-esima
    end
    
    [best_mesh] = Mesh_select(path_mesh, br_to_nod);   %%% funzione che permette di troncare il percorso delle maglie al primo nodo comune incontrato (procedendo a ritroso), senza arrivare sino al nodo di slack
    
    nod_idx_mesh_vect = nod_idx_mesh(j)*ones(size(best_mesh,1),1);
    best_mesh = [nod_idx_mesh_vect, best_mesh];        %%% si aggiunge il nodo finale di arrivo per la maglia, ovvero quello comune in cui convergono i due o più percorsi
    
    if size(best_mesh,2) > size(path_best_mesh,2)      %%% la parte di questo if consente di creare gli zeri da aggiungere per avere una matrice rettangolare quando vengono aggiunti più percorsi relativi alle maglie
       num_zeros_col = size(best_mesh,2) - size(path_best_mesh,2);
       num_zeros_rig = size(path_best_mesh,1);
       path_best_mesh = [path_best_mesh, zeros(num_zeros_rig, num_zeros_col)];
    end
    if size(path_best_mesh,2) > size(best_mesh,2)
        num_zeros_col = size(path_best_mesh,2) - size(best_mesh,2);
        num_zeros_rig = size(best_mesh,1);
        best_mesh = [best_mesh, zeros(num_zeros_rig, num_zeros_col)];
    end
    
    path_best_mesh = [path_best_mesh; best_mesh];
    max_dim = 0;
    for i=1:size(path_best_mesh,1)                     %%% la parte di questo for consente di eliminare tutte le colonne con elementi tutti nulli
        max_dim2 = find(path_best_mesh(i,:),1,'last');
        if max_dim2 > max_dim
            max_dim = max_dim2;
        end
    end
    path_best_mesh = path_best_mesh(:,1:max_dim);     %%% percorsi minimi che formano le maglie (saranno 2 per ciascuna maglia)
end

num_way = size(path_best_mesh,1);
rig = size(path_best_mesh,1);
col = size(path_best_mesh,2);
impedance_best_mesh = cell(rig, col-1);
branch_best_mesh = zeros(rig, col-1);

for i=1:num_way
    num_nod_way = find(path_best_mesh(i,:),1,'last');
    way = path_best_mesh(i,1:num_nod_way);
    br_cod = [];
    for j=1:length(way)-1
        idx_br = branch_end == way(1,j) & branch_first == way(1,j+1);
        if j == 1
            Zbr = Zbr_set(idx_br,1);
        else
            Zbr = [Zbr, Zbr_set(idx_br,1)];
        end
        br_cod = [br_cod, branch_cod(idx_br)];
    end
    impedance_best_mesh(i,1:length(Zbr)) = Zbr;      %%% impedenze associate ai rami corrispondenti ai percorsi minimi delle maglie (path_best_mesh)
    branch_best_mesh(i,1:length(br_cod)) = br_cod;   %%% codice dei rami corrispondenti ai percorsi minimi delle maglie (path_best_mesh)
end

num_rig = size(impedance_best_mesh,1);
num_rig = num_rig/2;
GzA = zeros(2*num_rig,6*num_branches); GzB = GzA; GzC = GzA;
GzcomplexA = zeros(num_rig, 3*num_branches); GzcomplexB = GzcomplexA; GzcomplexC = GzcomplexA;

for i=1:num_rig
    path1 = branch_best_mesh(2*i-1,:);       %%% percorso 1 relativo alla maglia i considerata
    path2 = branch_best_mesh(2*i,:);         %%% percorso 2 relativo alla maglia i considerata
    br1 = find(path1);
    br2 = find(path2);
    imp1 = impedance_best_mesh(2*i-1,:);
    imp2 = impedance_best_mesh(2*i,:);
    path1 = path1(br1);                      %%% percorso 1 senza elementi nulli
    path2 = path2(br2);                      %%% percorso 2 senza elementi nulli
    impedance1 = imp1(br1);                  %%% impedenze complesse associate al percorso 1
    impedance2 = imp2(br2);                  %%% impedenze complesse associate al percorso 2
    for j = 1:size(impedance1,2)
        k = path1(1,j);
        GzA(i, k) = real(impedance1{1,j}(1,1));
        GzA(i, k + num_branches) = real(impedance1{1,j}(1,2));
        GzA(i, k + 2*num_branches) = real(impedance1{1,j}(1,3));
        GzA(i, k + 3*num_branches) = - imag(impedance1{1,j}(1,1));
        GzA(i, k + 4*num_branches) = - imag(impedance1{1,j}(1,2));
        GzA(i, k + 5*num_branches) = - imag(impedance1{1,j}(1,3));
        GzB(i, k) = real(impedance1{1,j}(2,1));
        GzB(i, k + num_branches) = real(impedance1{1,j}(2,2));
        GzB(i, k + 2*num_branches) = real(impedance1{1,j}(2,3));
        GzB(i, k + 3*num_branches) = - imag(impedance1{1,j}(2,1));
        GzB(i, k + 4*num_branches) = - imag(impedance1{1,j}(2,2));
        GzB(i, k + 5*num_branches) = - imag(impedance1{1,j}(2,3));
        GzC(i, k) = real(impedance1{1,j}(3,1));
        GzC(i, k + num_branches) = real(impedance1{1,j}(3,2));
        GzC(i, k + 2*num_branches) = real(impedance1{1,j}(3,3));
        GzC(i, k + 3*num_branches) = - imag(impedance1{1,j}(3,1));
        GzC(i, k + 4*num_branches) = - imag(impedance1{1,j}(3,2));
        GzC(i, k + 5*num_branches) = - imag(impedance1{1,j}(3,3));
        
        GzA(i+num_rig, k) = imag(impedance1{1,j}(1,1));
        GzA(i+num_rig, k+num_branches) = imag(impedance1{1,j}(1,2));
        GzA(i+num_rig, k+2*num_branches) = imag(impedance1{1,j}(1,3));
        GzA(i+num_rig, k+3*num_branches) = real(impedance1{1,j}(1,1));
        GzA(i+num_rig, k+4*num_branches) = real(impedance1{1,j}(1,2));
        GzA(i+num_rig, k+5*num_branches) = real(impedance1{1,j}(1,3));
        GzB(i+num_rig, k) = imag(impedance1{1,j}(2,1));
        GzB(i+num_rig, k+num_branches) = imag(impedance1{1,j}(2,2));
        GzB(i+num_rig, k+2*num_branches) = imag(impedance1{1,j}(2,3));
        GzB(i+num_rig, k+3*num_branches) = real(impedance1{1,j}(2,1));
        GzB(i+num_rig, k+4*num_branches) = real(impedance1{1,j}(2,2));
        GzB(i+num_rig, k+5*num_branches) = real(impedance1{1,j}(2,3));
        GzC(i+num_rig, k) = imag(impedance1{1,j}(3,1));
        GzC(i+num_rig, k+num_branches) = imag(impedance1{1,j}(3,2));
        GzC(i+num_rig, k+2*num_branches) = imag(impedance1{1,j}(3,3));
        GzC(i+num_rig, k+3*num_branches) = real(impedance1{1,j}(3,1));
        GzC(i+num_rig, k+4*num_branches) = real(impedance1{1,j}(3,2));
        GzC(i+num_rig, k+5*num_branches) = real(impedance1{1,j}(3,3));
        
        GzcomplexA(i, k) = impedance1{1,j}(1,1);
        GzcomplexA(i, k+num_branches) = impedance1{1,j}(1,2);
        GzcomplexA(i, k+2*num_branches) = impedance1{1,j}(1,3);
        GzcomplexB(i, k) = impedance1{1,j}(2,1);
        GzcomplexB(i, k+num_branches) = impedance1{1,j}(2,2);
        GzcomplexB(i, k+2*num_branches) = impedance1{1,j}(2,3);
        GzcomplexC(i, k) = impedance1{1,j}(3,1);
        GzcomplexC(i, k+num_branches) = impedance1{1,j}(3,2);
        GzcomplexC(i, k+2*num_branches) = impedance1{1,j}(3,3); 
    end
    
    for j = 1:size(impedance2,2)
        k = path2(1,j);
        GzA(i, k) = - real(impedance2{1,j}(1,1));
        GzA(i, k + num_branches) = - real(impedance2{1,j}(1,2));
        GzA(i, k + 2*num_branches) = - real(impedance2{1,j}(1,3));
        GzA(i, k + 3*num_branches) = imag(impedance2{1,j}(1,1));
        GzA(i, k + 4*num_branches) = imag(impedance2{1,j}(1,2));
        GzA(i, k + 5*num_branches) = imag(impedance2{1,j}(1,3));
        GzB(i, k) = - real(impedance2{1,j}(2,1));
        GzB(i, k + num_branches) = - real(impedance2{1,j}(2,2));
        GzB(i, k + 2*num_branches) = - real(impedance2{1,j}(2,3));
        GzB(i, k + 3*num_branches) = imag(impedance2{1,j}(2,1));
        GzB(i, k + 4*num_branches) = imag(impedance2{1,j}(2,2));
        GzB(i, k + 5*num_branches) = imag(impedance2{1,j}(2,3));
        GzC(i, k) = - real(impedance2{1,j}(3,1));
        GzC(i, k + num_branches) = - real(impedance2{1,j}(3,2));
        GzC(i, k + 2*num_branches) = - real(impedance2{1,j}(3,3));
        GzC(i, k + 3*num_branches) = imag(impedance2{1,j}(3,1));
        GzC(i, k + 4*num_branches) = imag(impedance2{1,j}(3,2));
        GzC(i, k + 5*num_branches) = imag(impedance2{1,j}(3,3));
        
        GzA(i+num_rig, k) = - imag(impedance2{1,j}(1,1));
        GzA(i+num_rig, k+num_branches) = - imag(impedance2{1,j}(1,2));
        GzA(i+num_rig, k+2*num_branches) = - imag(impedance2{1,j}(1,3));
        GzA(i+num_rig, k+3*num_branches) = - real(impedance2{1,j}(1,1));
        GzA(i+num_rig, k+4*num_branches) = - real(impedance2{1,j}(1,2));
        GzA(i+num_rig, k+5*num_branches) = - real(impedance2{1,j}(1,3));
        GzB(i+num_rig, k) = - imag(impedance2{1,j}(2,1));
        GzB(i+num_rig, k+num_branches) = - imag(impedance2{1,j}(2,2));
        GzB(i+num_rig, k+2*num_branches) = - imag(impedance2{1,j}(2,3));
        GzB(i+num_rig, k+3*num_branches) = - real(impedance2{1,j}(2,1));
        GzB(i+num_rig, k+4*num_branches) = - real(impedance2{1,j}(2,2));
        GzB(i+num_rig, k+5*num_branches) = - real(impedance2{1,j}(2,3));
        GzC(i+num_rig, k) = - imag(impedance2{1,j}(3,1));
        GzC(i+num_rig, k+num_branches) = - imag(impedance2{1,j}(3,2));
        GzC(i+num_rig, k+2*num_branches) = - imag(impedance2{1,j}(3,3));
        GzC(i+num_rig, k+3*num_branches) = - real(impedance2{1,j}(3,1));
        GzC(i+num_rig, k+4*num_branches) = - real(impedance2{1,j}(3,2));
        GzC(i+num_rig, k+5*num_branches) = - real(impedance2{1,j}(3,3));
        
        GzcomplexA(i, k) = - impedance2{1,j}(1,1);
        GzcomplexA(i, k+num_branches) = - impedance2{1,j}(1,2);
        GzcomplexA(i, k+2*num_branches) = - impedance2{1,j}(1,3);
        GzcomplexB(i, k) = - impedance2{1,j}(2,1);
        GzcomplexB(i, k+num_branches) = - impedance2{1,j}(2,2);
        GzcomplexB(i, k+2*num_branches) = - impedance2{1,j}(2,3);
        GzcomplexC(i, k) = - impedance2{1,j}(3,1);
        GzcomplexC(i, k+num_branches) = - impedance2{1,j}(3,2);
        GzcomplexC(i, k+2*num_branches) = - impedance2{1,j}(3,3);        
    end
end

if size(GzA,1) == 0            %%% se non ci sono maglie si fissano dei vettori nulli per non avere problemi nella normale esecuzione dell'algoritmo
    GzA = [];
    GzcomplexA = [];
end
if size(GzB,1) == 0            %%% se non ci sono maglie si fissano dei vettori nulli per non avere problemi nella normale esecuzione dell'algoritmo
    GzB = [];
    GzcomplexB = [];
end
if size(GzC,1) == 0            %%% se non ci sono maglie si fissano dei vettori nulli per non avere problemi nella normale esecuzione dell'algoritmo
    GzC = [];
    GzcomplexC = [];
end

isA = nominal_linedata(:,13);
isB = nominal_linedata(:,14);
isC = nominal_linedata(:,15);

num_branchesA = sum(isA==1);
num_branchesB = sum(isB==1);
num_branchesC = sum(isC==1);

offsetB = num_branchesA;
offsetC = num_branchesB + offsetB;
offsetAx = num_branchesC + offsetC;
offsetBx = num_branchesA + offsetAx;
offsetCx = num_branchesB + offsetBx;

Gz = [GzA; GzB; GzC];
Gzcomplex = [GzcomplexA; GzcomplexB; GzcomplexC];
dim = size(Gz,1)/3;
dim2 = size(Gzcomplex,1)/3;
if dim >0
    isABC1 = [isA; isB; isC; isA; isB; isC]';
    isABC_logic1 = isABC1 == 1;
    Gz = Gz(:,isABC_logic1);

    isABC2 = [isA; isB; isC];
    isABC_logic2 = isABC2 == 1;
    Gzcomplex = Gzcomplex(:,isABC_logic2);
    Gzmod = abs(Gzcomplex);
    Gzangle = angle(Gzcomplex);

    GzA_state = [Gz(1:dim,1:offsetB),Gz(1:dim,offsetAx+1:offsetBx)];
    GzB_state = [Gz(dim+1:2*dim,offsetB+1:offsetC),Gz(dim+1:2*dim,offsetBx+1:offsetCx)];
    GzC_state = [Gz(2*dim+1:end,offsetC+1:offsetAx), Gz(2*dim+1:end,offsetCx+1:end)];
    Gz_state = cell(1,3);
    Gz_state{1,1} = GzA_state;
    Gz_state{1,2} = GzB_state;
    Gz_state{1,3} = GzC_state;
    
    GzmodA_state = Gzmod(1:dim2,1:offsetB);
    GzmodB_state = Gzmod(dim2+1:2*dim2,offsetB+1:offsetC);
    GzmodC_state = Gzmod(2*dim2+1:end,offsetC+1:end);
    GzangleA_state = Gzangle(1:dim2,1:offsetB);
    GzangleB_state = Gzangle(dim2+1:2*dim2,offsetB+1:offsetC);
    GzangleC_state = Gzangle(2*dim2+1:end,offsetC+1:end);
    
    Gzmod_state = cell(1,3);
    Gzangle_state = cell(1,3);
    Gzmod_state{1,1} = GzmodA_state;
    Gzmod_state{1,2} = GzmodB_state;
    Gzmod_state{1,3} = GzmodC_state;
    Gzangle_state{1,1} = GzangleA_state;
    Gzangle_state{1,2} = GzangleB_state;
    Gzangle_state{1,3} = GzangleC_state;
    
else
    Gz_state = cell(1,3);
    Gzmod_state = cell(1,3);
    Gzangle_state = cell(1,3);
    Gzmod = [];
    Gzangle = [];
end

Gzcomplex2 = Gzcomplex;
eq_curr_matrix2 = eq_curr_matrix(7:end,7:end);
dim_eq1 = size(eq_curr_matrix2,1);
dim_eq2 = size(eq_curr_matrix2,2);
eq_curr_matrix2 = eq_curr_matrix2(1:dim_eq1/2,1:dim_eq2/2);
if size(Gzcomplex,1) ~= 0;
    dim = size(Gzcomplex,1)/3;
    for i = 1:dim
        variab = branch_best_mesh((i*2)-1:i*2,1);
        idx = max(variab);
        index = idx - num_nodes + 1;    
        Gzcomplex2(index,:) = Gzcomplex(i,:);
        Gzcomplex2(index+dim,:) = Gzcomplex(i+dim,:);
        Gzcomplex2(index+2*dim,:) = Gzcomplex(i+2*dim,:);
    end
    Gzcomplex2 = Gzcomplex2*eq_curr_matrix2;
    Gzcomplex3 = Gzcomplex2;
    dim2 = size(Gzcomplex2,2);
    eq_matr = eye(dim2);
    dim3 = size(Gzcomplex,1);
    for i = 1:dim
        index = 3*dim + 1 -i;
        index2 = dim2 +1 -i;
        eq_matr(index2,:) = -Gzcomplex3(index,:)/Gzcomplex3(index,index2);
        eq_matr(index2,index2) = 0;
        eq_molt = zeros(dim2);
        eq_molt(index2,:) = eq_matr(index2,:);
        eq_matr = eq_matr + eq_matr*eq_molt;
        eq_matr(:,index2) = zeros(dim2,1);
        Gzcomplex3 = Gzcomplex3 + Gzcomplex3*eq_molt;
        Gzcomplex3(:,index2) = zeros(dim3,1);
    end
    for i = 1:dim
        index = 2*dim + 1 -i;
        index2 = dim2 +1 -i - num_branchesC_red;
        eq_matr(index2,:) = -Gzcomplex3(index,:)/Gzcomplex3(index,index2);
        eq_matr(index2,index2) = 0;
        eq_molt = zeros(dim2);
        eq_molt(index2,:) = eq_matr(index2,:);
        eq_matr = eq_matr + eq_matr*eq_molt;
        eq_matr(:,index2) = zeros(dim2,1);
        Gzcomplex3 = Gzcomplex3 + Gzcomplex3*eq_molt;
        Gzcomplex3(:,index2) = zeros(dim3,1);
    end
    for i = 1:dim
        index = dim + 1 -i;
        index2 = dim2 +1 -i - num_branchesC_red - num_branchesB_red;
        eq_matr(index2,:) = -Gzcomplex3(index,:)/Gzcomplex3(index,index2);
        eq_matr(index2,index2) = 0;
        eq_molt = zeros(dim2);
        eq_molt(index2,:) = eq_matr(index2,:);
        eq_matr = eq_matr + eq_matr*eq_molt;
        eq_matr(:,index2) = zeros(dim2,1);
        Gzcomplex3 = Gzcomplex3 + Gzcomplex3*eq_molt;
        Gzcomplex3(:,index2) = zeros(dim3,1);
    end
    eq_matr_r = real(eq_matr);
    eq_matr_x = imag(eq_matr);
    eq_tot = [eq_matr_r, -eq_matr_x; eq_matr_x, eq_matr_r];
    eq_matr2 = blkdiag(eye(6),eq_tot);
    eq_matr2(:,2*dim2-dim+7:2*dim2+6) = [];
    eq_matr2(:,2*dim2-dim+7-num_branchesC_red:2*dim2+6-num_branchesC_red) = [];
    eq_matr2(:,2*dim2-dim+7-num_branchesC_red-num_branchesB_red:2*dim2+6-num_branchesC_red-num_branchesB_red) = [];
    eq_matr2(:,dim2+7-dim:dim2+6) = [];
    eq_matr2(:,dim2+7-dim-num_branchesC_red:dim2+6-num_branchesC_red) = [];
    eq_matr2(:,dim2+7-dim-num_branchesC_red-num_branchesB_red:dim2+6-num_branchesC_red-num_branchesB_red) = [];
    eq_new = eq_curr_matrix*eq_matr2;
    eq_curr_matrix = eq_new;
    eq_curr_matrix = sparse(eq_curr_matrix);
end

