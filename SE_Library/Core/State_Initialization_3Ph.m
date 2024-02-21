function [V_abs, V_theta, I_abs, I_theta] = State_Initialization_3Ph(num_nodes, nominal_linedata, node_3ph_reverse, zdata, Zpaths, V_abs, V_theta, index, Gzcomplex)
    
    num_branches = size(nominal_linedata,1);
    
    branch_first_A = nominal_linedata(:,3);
    branch_first_B = nominal_linedata(:,4);
    branch_first_C = nominal_linedata(:,5);
    
    branch_end_A = nominal_linedata(:,6);
    branch_end_B = nominal_linedata(:,7);
    branch_end_C = nominal_linedata(:,8);
    
    num_branchesA = sum(branch_first_A > 0);
    num_branchesB = sum(branch_first_B > 0);
    num_branchesC = sum(branch_first_C > 0);
    
    A_branches_indices = find(branch_first_A > 0);
    B_branches_indices = find(branch_first_B > 0);
    C_branches_indices = find(branch_first_C > 0);
    
    A_branch_first = branch_first_A(A_branches_indices);
    B_branch_first = branch_first_B(B_branches_indices);
    C_branch_first = branch_first_C(C_branches_indices);
    
    A_branch_end = branch_end_A(A_branches_indices);
    B_branch_end = branch_end_B(B_branches_indices);
    C_branch_end = branch_end_C(C_branches_indices);
    
    num_A_nodes = max(branch_end_A);
    num_B_nodes = max(branch_end_B);
    num_C_nodes = max(branch_end_C);
   
    nodes_with_A = node_3ph_reverse(1:num_A_nodes,1);
    nodes_with_A = nodes_with_A(2:end); % escludo il primo che non riguardera le potenze
    nodes_with_B = node_3ph_reverse(1:num_B_nodes,2);
    nodes_with_B = nodes_with_B(2:end); % escludo il primo che non riguardera le potenze
    nodes_with_C = node_3ph_reverse(1:num_C_nodes,3);
    nodes_with_C = nodes_with_C(2:end); % escludo il primo che non riguardera le potenze
    
    nmeshA = size(Gzcomplex,1)/3; nmeshB = size(Gzcomplex,1)/3; nmeshC = size(Gzcomplex,1)/3;
    GzcomplexA = Gzcomplex(1:nmeshA,:); GzcomplexB = Gzcomplex(nmeshA+1:nmeshA+nmeshB,:); GzcomplexC = Gzcomplex(nmeshA+nmeshB+1:end,:);
    
    V = [V_abs(index);V_abs(index + num_nodes);V_abs(index + 2*num_nodes)];
    theta = [V_theta(index);V_theta(index+num_nodes);V_theta(index+2*num_nodes)];
    V_complex = V.*exp(1i*theta);

    SmatrixAA = zeros(num_branchesA, num_branchesA);
    SmatrixAB = zeros(num_branchesA, num_branchesB);
    SmatrixAC = zeros(num_branchesA, num_branchesC);
    SmatrixBA = zeros(num_branchesB, num_branchesA);
    SmatrixBB = zeros(num_branchesB, num_branchesB);
    SmatrixBC = zeros(num_branchesB, num_branchesC);
    SmatrixCA = zeros(num_branchesC, num_branchesA);
    SmatrixCB = zeros(num_branchesC, num_branchesB);
    SmatrixCC = zeros(num_branchesC, num_branchesC);
    
    for i=1:num_branchesA
        fbA = A_branch_first(i);
        tbA = A_branch_end(i);
        if (fbA ~= 1)
            SmatrixAA(fbA-1,i) = V_complex(1);
        end
        SmatrixAA(tbA-1,i) = -V_complex(1);
    end
    for i=1:num_branchesB
        fbB = B_branch_first(i);
        tbB = B_branch_end(i);
        if (fbB ~= 1)
            SmatrixBB(fbB-1,i) = V_complex(2);
        end
        SmatrixBB(tbB-1,i) = -V_complex(2);
    end
    for i=1:num_branchesC
        fbC = C_branch_first(i);
        tbC = C_branch_end(i);
        if (fbC ~= 1)
            SmatrixCC(fbC-1,i) = V_complex(3);
        end
        SmatrixCC(tbC-1,i) = -V_complex(3);
    end
    
    
    SmatrixA = [SmatrixAA, SmatrixAB, SmatrixAC];
    SmatrixB = [SmatrixBA, SmatrixBB, SmatrixBC];
    SmatrixC = [SmatrixCA, SmatrixCB, SmatrixCC];
    
    if size(Gzcomplex,1) > 0
        SmatrixA(end-nmeshA+1:end,:) = conj(GzcomplexA);
        SmatrixB(end-nmeshB+1:end,:) = conj(GzcomplexB);
        SmatrixC(end-nmeshC+1:end,:) = conj(GzcomplexC);
    end
    
    type = zdata(:,1);
    pii = (type == 2);
    qii = (type == 3);
    
    SA = complex(zdata(pii,2),zdata(qii,2));
    SB = complex(zdata(pii,3),zdata(qii,3));
    SC = complex(zdata(pii,4),zdata(qii,4));
    SA = SA(nodes_with_A-1);
    SB = SB(nodes_with_B-1);
    SC = SC(nodes_with_C-1);
    SA = [SA; zeros(nmeshA,1)];
    SB = [SB; zeros(nmeshB,1)];
    SC = [SC; zeros(nmeshC,1)];
    
    Smatrix = [SmatrixA; SmatrixB; SmatrixC];
    S = [SA; SB; SC];
    Ibr = conj(Smatrix) \ conj(S);
    index = [A_branches_indices; B_branches_indices + num_branches; C_branches_indices + 2*num_branches];
    I = zeros(3*num_branches,1);
    I(index) = Ibr;
    
    I_abs = abs(I);
    I_theta = angle(I);
    
    V = kron(V_complex,ones(num_nodes,1)) - Zpaths * I;
    
    V_abs = abs(V);
    V_theta = angle(V);

end