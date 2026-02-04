README

T3ST requires the following folder structure:

Parent folder
	- config
        - data
	    - Sim_01
		- Sim_02
		.....
	- G_EDSQK (potentially optional)
        - Mathematica (optional)
	- Sims       
	- Source
        T3ST_GUI.py  (GUI)
  	 
	 
	Order of steps: 
	    1. Terminal -> Parent folder -> python3 T3ST_GUI.py  (double click/enter works under Windows)
		2. The GUI Opens 
		3. Choose the number of the Simulation; 
		    Choose the scenario
			Modify manually parameters if needed
			Choose the parameter to be changed
			How many discrete values (no of points)
			What interval (init; final)
			Export parameters
		4. Goto: Parent folder/Source
		        source Instructiuni_0 + source Instructiuni_1
		        make (alternative but there are some unclear issues)
			./a.out
			Input the no of the simulation (noofsim)
		5. Process data from Parent folder/data/Sim_noofsim/Run_000??	(with Mathematica)
	    
			
