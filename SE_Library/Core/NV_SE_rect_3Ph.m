function [Vr, Vx, Ir, Ix, num_iter] = NV_SE_rect_3Ph(num_nodes, num_branches, nominal_linedata, node_3ph, node_3ph_reverse, V_abs, V_theta, zdata, Y_state, Y_tot, Znod, Gzcomplex)

type = zdata(:,1);
from = zdata(:,6);

Vmag_measi = (type == 1);  %%% Voltage amplitude measurements initialization 
Vmag_meas = find(Vmag_measi);

if length(Vmag_meas) >= 1
    idx = Vmag_meas(1);
    moduloA = zdata(idx,2);
    moduloB = zdata(idx,3);
    moduloC = zdata(idx,4);
    trad_code = 1;
    if (moduloA == 0)
        if (moduloB == 0)
            moduloA = moduloC;
            moduloB = moduloC;
        else
            moduloA = moduloB;
            if (moduloC == 0)
                moduloC = moduloB;
            end
        end
    else
        if (moduloB == 0)
            moduloB = moduloA;
        end
        if (moduloC == 0)
            moduloC = moduloA;
        end
    end
    modulo = [moduloA; moduloB; moduloC];
    fase = [0,0,0];
    fase = kron(fase',ones(num_nodes,1));
    index = from(idx);
end

Vsync_measi = (type == 7);                        
Vsync_meas = find(Vsync_measi);
if length(Vsync_meas) >= 1
    disp([mfilename, ': synchronized measurements management not available']);
end

Zpaths = zeros(3*num_nodes,3*(num_branches));
indicesForAllButFirstNode = 1 : 3*num_nodes; 
indicesForAllButFirstNode([1, num_nodes+1, 2*num_nodes + 1]) = [];
Zpaths(indicesForAllButFirstNode,1:end) = Znod;

ZrotA = ones(num_nodes,1)*Zpaths(index,:);
ZrotB = ones(num_nodes,1)*Zpaths(index + num_nodes,:);
ZrotC = ones(num_nodes,1)*Zpaths(index + 2*+ num_nodes,:);
Zpaths_ini = Zpaths - [ZrotA; ZrotB; ZrotC];
V_abs_in = kron(modulo,ones(num_nodes,1));

[V_abs, V_theta, I_abs, I_theta] = State_Initialization_3Ph(num_nodes, nominal_linedata, node_3ph_reverse, zdata, Zpaths_ini, V_abs_in, V_theta+fase, index, Gzcomplex);


V_theta(1:num_nodes) = V_theta(1:num_nodes) - ones(num_nodes,1)*V_theta(1);
V_theta(num_nodes+1:2*num_nodes) = V_theta(num_nodes+1:2*num_nodes) - ones(num_nodes,1)*V_theta(num_nodes+1) + 2*pi/3;
V_theta(2*num_nodes+1:end) = V_theta(2*num_nodes+1:end) - ones(num_nodes,1)*V_theta(2*num_nodes+1) - 2*pi/3;
[Vr, Vx, Ir, Ix, num_iter] = NV_DSSE_rectangular_traditional_3Ph(num_nodes, nominal_linedata, node_3ph, node_3ph_reverse, V_abs, V_theta, I_abs, I_theta, zdata, Y_state, Y_tot);

end

