De-immunization of Factor VIII using the evolutionary Hamiltonian
========

Supplementing code of article: De-immunization of Factor VIII using the evolutionary Hamiltonian

Authors: Benjamin Schubert 
Date: Mai 2015  
Version: 1.0  
License: The software is released under a three-clause BSD license


Introduction:
-------------
The provided software enables the solving of bi-objective mixed integer problems (as long as only one objective is mixed integer) such as those
problems arising in the corresponding paper. The solver is a parallel two-phase extension of the rectangle-splitting approach proposed by Boland et al. 2015. Detail description of the algorithm can be found in the supplementary material. The approach is used to solve the proposed de-immunization problem, but can also be used as generic BOMIP-Solver (as long as only one objective is mixed integer).


Requirement:
-------------
The solver uses the following software and libraries:  

1. Python 2.7
2. Numpy 1.9.1
3. Cplex (+Python API) 12.5
4. Polygon 2.0.7

Please make sure you have installed said software/libraries
and their dependencies.


Installation:
-------------
Please install all required libraries. The Python API of CPLEX should be installed as well. CPLEX is free for academic use. For more 
details see IBMs Academic Initiative (http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative).

Structure:
-------------
- AMPL/continuous_immuno_hemiltonian_en.mod (AMPL formulation of the BOMIP de-immunization problem - Energy objective)
- AMPL/continuous_immuno_hemiltonian_imm.mod (AMPL formulation of the BOMIP de-immunization problem - Immunogenicity objective)
- src/deimmunization_manager.py 	  (Main Algorithm)
- src/deimmunization_worker.py 	    (Core solving step)
- src/utility/Solution.py 			 	  (Object storing the BOMIP solution, warmstart-configuration and variables of interrest)
- src/utility/Hypervolume.py 				(Functionality to approximate the quality of the current solution via Hypervolume calculation)
- src/utility/ParetoFilter.py 			(Functionality to filter duplicat and non-pareto optimal solutions)
- src/examples/                     (Instances solved in the main manuscript in cplex-lp format; solutions as cPickle python objects)


Usage:
-------------
The solver can be used in distributed systems. For that, first start the manager script, which implements the algorithm and handles the work distribution.
Make sure that the chosen IP and Port combination is reachable from all other used systems.
```
usage: deimmunization_manager.py [-h] [--grid GRID] --port PORT
                                 [--approximate APPROXIMATE] [--key KEY]
                                 --output OUTPUT [--resolve RESOLVE]

Rectangle Manager implementation

Arguments:
  -h, --help            show this help message and exit
  --grid GRID, -g GRID  Number of Epsilon grid points (default:3)
  --port PORT, -p PORT  Port
  --approximate APPROXIMATE, -a APPROXIMATE
                        Bound on approximation (default:0.009)
  --key KEY, -k KEY     Authentication key (default: rectangle)
  --output OUTPUT, -o OUTPUT
                        Solution output as pickel
  --resolve RESOLVE, -r RESOLVE
                        Reinitialize with partial solution
```

After initializing the manager process, one can start multiple worker processes on the same or distributed machines. The worker process connects to the manager process via TCP/IP at the specified port (must be the same as the one the manager listens at) and obtains single work packages from the manager process. Using the flag -r (--resolve) one can specify a intermediate solution and refine or restart the solving process from there. The specified intermediate solution must be a pickled list of Solution objects.

```
usage: deimmunization_worker.py [-h] --input INPUT1 INPUT2 --masterip MASTERIP
                                --port PORT --authkey AUTHKEY --threads
                                THREADS

Rectangle Worker Grid implementation

Arguments:
  -h, --help            show this help message and exit
  --input INPUT1 INPUT2, -i INPUT INPUT
                        model files
  --masterip MASTERIP, -m MASTERIP
                        The IP of the master node
  --port PORT, -p PORT  port to connect
  --authkey AUTHKEY, -a AUTHKEY
                        authentication key
  --threads THREADS, -t THREADS
                        nof of core
```

INPUT1 and INPUT2 are single-objective files (in CPLEX compatible formats) of the bi-objective problems, each containing one of the objective function and all constraints of the bi-objective problem. 

Examples:
-------------
Manager:
```
python deimmunization_manager.py -p 6882 -g 3 -a 0.009 -k rectangle -o ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_output.pcl
```
Worker:
```
python deimmunization_worker.py -i ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_pssm05_f01_f_he_experiment_local_imm.lp ./example/FA8_HUMAN_hmmer_plm_n5_m50_f70_t01_g_r2188-2345_e20_k1_pssm05_f01_f_he_experiment_local_en.lp -m 127.0.0.1 -p 6882 -a rectangle -t 4
```

Contacts:
-------------
Benjamin Schubert  
schubert@informatik.uni-tuebingen.de  
University of Tübingen, Applied Bioinformatics,  
Center for Bioinformatics, Quantitative Biology Center,  
and Dept. of Computer Science,  
Sand 14, 72076 Tübingen, Germany
