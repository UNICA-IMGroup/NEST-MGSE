# NEST-MGSE
The software realizes a framework for A framework for microgid state estimation (MGSE). 
State estimation is a fundamental tool to augment the situational awareness of the system enabling the development of strategic applications.

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
			MGSE_Test_MC
			MGSE_Test_MC Test_MC_config_3
	
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
			MGSE_Test_MC
			MGSE_Test_MC Test_MC_config_3
			

	
	
