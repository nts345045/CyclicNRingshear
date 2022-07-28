# CyclicNRingshear

core/ - core processing routines
	STEP0_Conversion/ - for raw data conversion and instrument corrections 
		RS_dataConversion_OSC_N.m - MATLAB script for converting transducer voltage measures into calibrated physical units 
		merge_lvdt.py - concatenates individual LVDT files and applies corrections for instrument position shifts 
		unify_timestamp_indexing.py - Convert all data into local datetime indexing (UTC - 5 hours) 
	STEP1_Cleaning/ -  for cleaning time-series 
		smooth_despike_and_segment.py - Apply post-processing to signals and segment into experimental/steady-state subsets 
										also generates Fig. S1
		process_cavity_picks.py - Process detachment/crest/reattachment point geometries picked from images to create  
								  timeseries of contact geometries
			process_LVDT.py - conduct melt-rate corrections on experimental and steady-state data (some results written to command line)
	STEP2_Interpretation/ - for data interpretation and figure rendering
		CNSB_Fig2_Forcing.py - Create figure from effective pressure data that marks steady-state and experimental periods & cycle numbers 
		CNSB_Figs3and4_T24_Timeseries.py - Create figures of timeseries data from experiment T24 and end of steady-state period
		CNSB_Fig6_T6_Timeseries.py - Create figure of timeseries data from experiment T6 
		CNSB_Fig7_T96_Timeseries.py - Create figure of timeseries data from experiment T96 
		CNSB_Fig8_Hysteresis.py - Create figure of hysteresis loops for all three experiments 
		CNSB_Fig9_NormModel.py - Create figure comparing data from experiments in T24 and T6 to steady-state theory 
util/ - Supporting scripts for core processing
	datetimeindex.py - methods for working with pandas.DatetimeIndex objects and associated data vectors	