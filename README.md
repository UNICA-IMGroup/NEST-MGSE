# NEST-MGSE
The software realizes a for microgid state estimation (MGSE). 
State estimation is a fundamental tool to augment the situational awareness of the system enabling the development of strategic applications.
State Estimation is based on mathematical relations between system state variables and available measurements that can be generally represented by means of the so-called measurement equation:
	z=h(x)+e	(1)
where z=[z_1,z_2,…,z_M ]^T is the vector of the M available measurements, h is the mathematical relation representing the measurement functions (depending on the specific measurement type, 
for example, power flow or injection, voltage or current measurements, etc.) that involve network topology and parameters, x=[x_1,x_2,…,x_N ]^T is the vector of the N state variables and e is the measurement error vector. 
In the case of poorly monitored networks, as in the cases addressed in this study, it is possible to exploit to so-called pseudo-measurements, i.e. a priori knowledge of the powers absorbed or injected into the various nodes of the network.
Stranting from the Network under monitoring description, the measurement system description, and the measurements (psedudo-measurements and real time measurements), the software provide the estimates of all the node Voltages and of all the branch currents, in module and phase angle, along with their extended uncertanty.  

Software Organization
	SE_Library - Folder
		It is the folder that contains the library.
		Library's interface routine:
			-Network_parameters_to_SE.m 
			-MeasurementConfiguration_to_SE.m
			-State_Estimation.m
		Core - Folder	
		Networks - Folder
			It contains Networks definition files
		Measurements - Folder
			-Measurement_placement.m it defines the measurement placement configurations
			-Measurement_uncertainty.m it defines the measurement uncertainty configurations
	SE_Tests - Folder
		It is the folder that performs the tests.
			-MGSE_Test.m is the main of the test, use:
				MGSE_Test
				MGSE_Test Test_MC_config_3
		
		SE_Utilities
			It contains utility files used in the tests, examples:
			-Newton_Raphson_Power_Flow_to_SE.m implements the power flow Newton-Raphson method to obtain reference values
			-AddMeasurementErrors_SE.m simulate measurements
			-PlotNetwork_to_SE.m plot the graph of the network
			-plot_MGSE(V, I, Network_param.topology, general_title); plot the estimate values with their confidence interval

How to use
SE_Tests - Folder
	It is the folder that performs the tests.
		-MGSE_Test.m is the main of the test, use:
			MGSE_Test
			MGSE_Test Test_MC_config_3
			

	
	
