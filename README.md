# **EarlyWorld** 
EarlyWorld is a computational environment that simulates, in a stochastic and discrete-time manner, the non-enzymatic polymerization and template-dependent replication of RNA on a clay-water interface. 
The system’s performance depends on two main parameters that represent the strength of two types of non-covalent interactions in a particular environment: 
* Clay - nucleotide (nt) interaction: Represented by parameter α.
* nt - nt interaction: Represented by parameter β.

EarlyWorld simulates two different compartments of the clay-water interface using independent simulations:
* Compartment I or Comp.I: Only single-stranded (ss) RNA molecules can form thanks to clay-nt interactions by random polymerization. To run a simulation of Comp.I use the files located in the folder EarlyWorld/CompI/.
* Compartment II or Comp.II: Both ssRNA but also double-stranded (ds) RNA molecules can form thanks to clay-nt and nt-nt interactions, respectively. Here random polymerization occurs, as in Comp.I, but the replication of an original molecule due to template-dependent polymerization also takes place. To run a simulation of Comp.II use the files located in the folder EarlyWorld/CompII/.

## Compartment I simulations 

### Main program (main.m)
The simulation starts with *N* single nucleotides floating in the aqueous phase and the clay level totally empty and reproduces the interaction of the clay and the nts in a particular environment, polymerization only taking place at the clay level. Depending on this environment, the process will give rise to the formation of oligomers that can remain adsorbed to the clay or dettach and return to the pool.

* Input:
  * **_A_**: String containing the different types of nts of the alphabet to use. E.g. A='augc'.
  * **_D_**: Matrix with the nucleobase complementarities allowed. E.g. D=\['au'; 'gc'\] means 'a' interacts with 'u' and 'g' interacts with 'c'.
  * **_N_**: Array of size 1-by-size(_A_), that contains the number of nucleotides of each type with which the simulation starts.
  * **_positions_**: Number of simulated clay positions.
  * **_t_**: Number of simulation time steps.
  * **_alpha_**: Array of size 1-by-_t_ that stores the value of parameter α to be used at each time step.
  * **_beta_**: Array of size 1-by-_t_ that stores the value of parameter β to be used at each time step. Notice: it is going to be an array of zeros for every simulation of Comp.I.

* Output (summary of the variables important for Comp.I simulations):
  * **_Activity_**: Cell that stores the total number of molecules that have been adsorbed and dettached from the clay over the course of the simulation, as well as this value split by polymer size.
  * **_Ct_**: String matrix of size _t_-by-_positions_ that stores the state of each clay position at each time step; i.e., wether it is empty or if it is occupied by a nt (specifying the nt type as in _A_).
  * **_Ot\_C_**: Matrix of size _t_-by-_positions_ that stores the state of each clay position at each time step: 0 if it is empty and 1 otherwise.
  * **_Lt_**: String matrix with _t_ rows that stores the state of the pool or aqueous phase at each time step; i.e., it contains all the molecules that are free in solution at each time step.
 
  ### Example
  To run an example simulation of Comp.I, use the EarlyWorld_example_simulation.m file as follows:
   ``` matlab
  EarlyWorld_example_simulation("Example","0.8","[0.5,2500]","1","100","1")
  ```
  Running the above code you obtain 100 stochastic realizations of an oscillating environment in which **_alpha_** follows $\alpha=0.8 + 0.5\sin{(2\pi t/2500)}$. The results are stored as 100 matlab files in the folder /Data/Example/.

## Compartment II simulations 

### Main program (main.m)
The simulation starts with *N* single nucleotides floating in the aqueous phase and one original ssRNA molecule with a random sequence placed at a random clay position and reproduces the interaction of the nts with the clay and with other nts (the latter via nucleobase complementarity) in a particular environment. Three levels are relevant for the evolution of the system: the clay or “level 0”; the complementary level or “level 1”, where one nt can only be placed if it is complementary to other nt included in a ssRNA oligomer which is placed in level 0; and the aqueous phase, pool or “level 2”. Depending on the environment, the process will give rise to the formation of copies of the original O sequence: both complementary C strands and replicate R strands (i.e., complementary to C, thus recovering the sequence of the original ssRNA molecule) form.

* Input:
  * **_A_**: String containing the different types of nts of the alphabet to use. E.g. A='augc'.
  * **_D_**: Matrix with the nucleobase complementarities allowed. E.g. D=\['au'; 'gc'\] means 'a' interacts with 'u' and 'g' interacts with 'c'.
  * **_N_**: Array of size 1-by-size(_A_), that contains the number of nucleotides of each type with which the simulation starts.
  * **_positions_**: Number of simulated clay positions.
  * **_t_**: Number of simulation time steps.
  * **_alpha_**: Array of size 1-by-_t_ that stores the value of parameter α to be used at each time step.
  * **_beta_**: Array of size 1-by-_t_ that stores the value of parameter β to be used at each time step.
  * **_polymer\_size_**: Length of the original ssRNA molecule.

* Output:
  * **_Activity_**: Cell that stores the total number of molecules that have been adsorbed or hybridized and dettached or denatured from level 0, level 1 and level 2 over the course of the simulation, as well as this value split by polymer size.
  * **_Ct_**: String matrix of size _t_-by-_positions_ that stores the state of each clay position at each time step; i.e., wether it is empty or if it is occupied by a nt (specifying the nt type as in _A_).
  * **_Ot\_C_**: Matrix of size _t_-by-_positions_ that stores the state of each clay position at each time step: 0 if it is empty and 1 otherwise.
  * **_Lt_**: String matrix with _t_ rows that stores the state of the pool or aqueous phase at each time step; i.e., it contains all the molecules that are free in solution at each time step.
  * **_Ut_**: String matrix of size _t_-by-_positions_ that stores the state of each complementary level position at each time step; i.e., wether it is empty or if it is occupied by a nt (specifying the nt type as in _A_).
  * **_Ot\_U_**: Matrix of size _t_-by-_positions_ that stores the state of each complementary level position at each time step: 0 if it is empty and 1 otherwise.
  * **_C\_CELL_**: Cell that stores all the C strands that have been formed and released back to the pool along the entire simulation and their characteristics.
  * **_R\_CELL_**: Cell that stores all the R strands that have been formed and released back to the pool along the entire simulation and their characteristics.

 * Other output variables:
In order to keep track of the identity of the oligomers formed in the simulation (i.e., wheter they are products of template-dependent polymerization or if they are random newly-generated molecules), other variables are created in the simulation and their values are stored. Those variables are: **_Lt\_OCR_**, **_Lt\_bl_**, **_Lt\_E_**, **_Ct\_OCR_**, **_Ct\_bl_**, **_Ct\_E_**, **_Ut\_OCR_** and **_Ut\_E_**. Their meaning is explained in the comments of the main.m function.

 ### Example
  To run an example simulation of Comp.II, use the EarlyWorld_example_simulation.m file as follows:
   ``` matlab
  EarlyWorld_example_simulation("Example","0.8","6","[0.5,5.9]","[2500,2500]","20","1","100","1","20000")
  ```
  Running the above code you obtain 100 stochastic realizations of an oscillating environment in which the original ssRNA molecule length is **20 nts**, **_alpha_** follows $\alpha=0.8 + 0.5\sin{(2\pi t/2500)}$ and **_beta_** follows $\beta=6 + 5.9\sin{(2\pi t/2500)}$. The results are stored as 100 matlab files in the folder /Data/Example/.

## Alphabet simulations 

These simulations proceed exactly in the same way as those of Comp.II, the only difference being the denaturation probability (rupture_alphabets.m), that in this case depends directly on the number of hydrogen bonds of each particular nt-nt interaction, instead of depending on the length of the molecule, as it was the case of Comp.II simulations.

 ### Example
  To run an example of alphabet simulations, use the EarlyWorld_example_alphabet.m file as follows:
   ``` matlab
  EarlyWorld_example_alphabet("Example","0.8","6","[0.5,5.9]","[2500,2500]","20","1","100","1","20000")
  ```
  Running the above code you obtain 100 stochastic realizations of an oscillating environment in which the original ssRNA molecule length is **20 nts**, the RNA **alphabet** is a two-letter one with G-C (three hydrogen bonds) nucleobase interactions, **_alpha_** follows $\alpha=0.8 + 0.5\sin{(2\pi t/2500)}$ and **_beta_** follows $\beta=6 + 5.9\sin{(2\pi t/2500)}$. The results are stored as 100 matlab files in the folder /Data/Example/.

## Hydrolysis simulations 

These simulations include a probability of hydrolysis per susceptible covalent bond ($P_{hyd}$), with the porpose of studying how it may affect replication. 

### Main program (main_hydrolysis.m)

The new variables that need to be incorporated into the main_hydrolysis function are the ones related to the calculation of hydrolysis probabilities. While the new output variables reflect the molecules hydrolyzed and the number of susceptible and broken covalent bonds.

* New input variable for the case in which $P_{hyd}$ is constant:
  * **_k_hyd_**: Value of the hydrolysis probability to be used in the simulation.
  
* New input variable for the case in which P_{hyd} oscillates:
  * **_k_hyd_**: Array of size 1-by-_t_ that stores the value of the hydrolysis probability to be used at each time step.
 
* New output variables (equal for both cases of constant and oscillatory hydrolysis probabilities):
  * **_HYDROLYZED_**: Cell that stores all the molecules that have been broken along the simulation, their labels, and the new molecules formed as a consequence of their hydrolysis.
  * **_N_pbonds_**: Array of size 1-by-_t_ that stores the number of susceptible covalent bonds at eacth time step.
  * **_H_pbonds_**: Array of size 1-by-_t_ that stores the number of hydrolyzed covalent bonds at eacth time step.


 ### Examples
  To run an example in which the probability of hydrolysis remains constant over time, use the EarlyWorld_hydrolysis.m file in the folder EarlyWorld/CompII/Hydrolysis/Constant/ as follows:
   ``` matlab
  EarlyWorld_hydrolysis("Example","0.8","6","[0.5,5.9]","[2500,2500]","20","1","100","1","20000","[1e-6]")
  ```
  Running the above code you obtain 100 stochastic realizations of an environment in which the original ssRNA molecule length is **20 nts**, **_alpha_** follows $\alpha=0.8 + 0.5\sin{(2\pi t/2500)}$, **_beta_** follows $\beta=6 + 5.9\sin{(2\pi t/2500)}$ and $P_{hyd}=10^{-6}$. The results are stored as 100 matlab files in the folder /Data/Example/.
  
To run an example in which the probability of hydrolysis oscillates over time, use the EarlyWorld_hydrolysis.m file in the folder EarlyWorld/CompII/Hydrolysis/Oscillatory/ as follows:
   ``` matlab
  EarlyWorld_hydrolysis("Example","0.8","6","[0.5,5.9]","[2500,2500]","20","1","100","1","20000","[-7]", "[-6]", "0")
  ```
  Running the above code you obtain 100 stochastic realizations of an oscillating environment in which the original ssRNA molecule length is **20 nts**, **_alpha_** follows $\alpha=0.8 + 0.5\sin{(2\pi t/2500)}$, **_beta_** follows $\beta=6 + 5.9\sin{(2\pi t/2500)}$ and $P_{hyd}$ fluctuates as follows: $\log P_{hyd}=\frac{1}{2}\log(10^{-7}\cdot 10^{-6})+\frac{1}{2}\log(10^{-7} / 10^{-6})\sin{(2\pi t / 2500)}$.
