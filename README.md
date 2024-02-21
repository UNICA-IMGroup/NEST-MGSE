# NEST-MGSE
The software realizes a framework for the state estimation in a generic distribution network as the first step to develop a state estimation for microgrids (MGSE). 
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
		-SE_Test_MC.m is the main of the test, use:
			SE_Test_MC
			SE_TEst_MC Test_MC_config_3
	
	SE_Utilities
		It contains utility files used in the tests, examples:
		-Newton_Raphson_Power_Flow_to_SE.m implements the power flow Newton-Raphson method to obtain reference values
		-AddMeasurementErrors_SE.m simulate measurements
		-PlotNetwork_to_SE.m plot the graph of the network
		-plotElapsedTime.m
		-plotVoltages_3Ph.m
		-plotCurrents_3Ph.m

How to use
SE_Tests - Folder
	It is the folder that performs the tests.
		-SE_Test_MC.m is the main of the test, use:
			SE_Test_MC
			SE_TEst_MC Test_MC_config_3
			

	
	