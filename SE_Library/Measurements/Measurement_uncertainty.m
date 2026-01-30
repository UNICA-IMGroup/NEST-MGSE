function meas_unc = Measurement_uncertainty(unc_config)
% input
% -unc_congig: is the code name 
% output
% -meas_unc: struct containg the measurements uncertainty
if ~exist('unc_config', 'var') 
  unc_config = 'default';  
end
switch unc_config
    case 'default'
        %Pseudo measurements
        meas_unc.unc_percent_PQ_pseudo = 50;
        
        % Standard measurement
        meas_unc.unc_percent_PQ_inj = 3;   
        meas_unc.unc_percent_PQ_br = 3;
        meas_unc.unc_percent_Vmag_nod = 1;
        meas_unc.unc_percent_Imag_br = 3;
  
    case 'nominal'
        %Pseudo measurements
        meas_unc.unc_percent_PQ_pseudo = 50;
        
        % Standard measurement
        meas_unc.unc_percent_PQ_inj = 1;   
        meas_unc.unc_percent_PQ_br = 1;
        meas_unc.unc_percent_Vmag_nod = 1;
        meas_unc.unc_percent_Imag_br = 1;

    case 'default + forecasting 10'
        %Pseudo measurements
        meas_unc.unc_percent_PQ_pseudo = 10;
        
        % Standard measurement
        meas_unc.unc_percent_PQ_inj = 3;   
        meas_unc.unc_percent_PQ_br = 3;
        meas_unc.unc_percent_Vmag_nod = 1;
        meas_unc.unc_percent_Imag_br = 3;
    
	case 'forecasting + 25'
        %Pseudo measurements
        meas_unc.unc_percent_PQ_pseudo = 25;
        
        % Standard measurement
        meas_unc.unc_percent_PQ_inj = 3;   
        meas_unc.unc_percent_PQ_br = 3;
        meas_unc.unc_percent_Vmag_nod = 1;
        meas_unc.unc_percent_Imag_br = 3;
        
    otherwise
        disp([mfilename, ': warning, unc_code not defined. default unc_code configuration is used']); 
        
        %Pseudo measurements
        meas_unc.unc_percent_PQ_pseudo = 50;
        
        % Standard measurement
        meas_unc.unc_percent_PQ_inj = 3;   
        meas_unc.unc_percent_PQ_br = 3;
        meas_unc.unc_percent_Vmag_nod = 1;
        meas_unc.unc_percent_Imag_br = 3;

end