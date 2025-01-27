function meas_indices = Measurement_placement(placement_config)
% input
% -placement_config is a code to define a specific placement. e.g 
%  'default' is used for the measurement configuration having only the
%  first node and the first branch monitored.
% output
% meas_indices, a struct of vectors that contains the indices of nodes and
% branches where the specific measurement is available
if ~exist('placement_config', 'var')
    placement_config = 'basic placement';
end

switch placement_config
    case 'default'    
        %Pseudo measurements
        %meas_indices.PQ_pseudo = [2 : num_nodes]';
        %Standard measurements
        meas_indices.PQ_br_idx = []';                                                              
        meas_indices.PQ_inj_idx = []';     
        meas_indices.Vmag_nod_idx = [1]'; 
        meas_indices.Imag_br_idx = [1]';
        %Synchonized measurements
%         meas_indices.Isync_magphase_br_idx = []';
%         meas_indices.Vsync_magphase_nod_idx = []';
    case '3 nodes placement'    
        %Pseudo measurements
        %meas_indices.PQ_pseudo = [2 : num_nodes]';
        %Standard measurements
        meas_indices.PQ_br_idx = []';                                                              
        meas_indices.PQ_inj_idx = []';     
        meas_indices.Vmag_nod_idx = [1, 7, 9]'; 
        meas_indices.Imag_br_idx = [1]';
        %Synchonized measurements
%         meas_indices.Isync_magphase_br_idx = []';
%         meas_indices.Vsync_magphase_nod_idx = []';
	case '4 nodes placement'    
			%Pseudo measurements
			%meas_indices.PQ_pseudo = [2 : num_nodes]';
			%Standard measurements
			meas_indices.PQ_br_idx = []';                                                              
			meas_indices.PQ_inj_idx = [7, 10, 13]';     
			meas_indices.Vmag_nod_idx = [1]'; 
			meas_indices.Imag_br_idx = [1]';
			%Synchonized measurements
	%         meas_indices.Isync_magphase_br_idx = []';
	%         meas_indices.Vsync_magphase_nod_idx = []';
    otherwise
        disp([mfilename, ': warning, placement configuration ', placement_config, ' not defined. basic placement configuration is used']);
        %Pseudo measurements
        %meas_indices.PQ_pseudo = [2 : num_nodes]';
        %Standard measurements
        meas_indices.PQ_br_idx = []';                                                               
        meas_indices.PQ_inj_idx = []';     
        meas_indices.Vmag_nod_idx = [1]'; 
        meas_indices.Imag_br_idx = [1]';
        %Synchonized measurements
%         meas_indices.Isync_magphase_br_idx = []';
%         meas_indices.Vsync_magphase_nod_idx = []';
end